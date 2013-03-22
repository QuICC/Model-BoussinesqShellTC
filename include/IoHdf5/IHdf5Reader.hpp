/** \file IHdf5Reader.hpp
 *  \brief Interface to a general HDF5 file reader
 */

#ifndef IHDF5READER_HPP
#define IHDF5READER_HPP

// Configuration includes
//
#include "Framework/FrameworkMacro.h"
#include "SmartPointers/SharedPtrMacro.h"

// System includes
//
#include <string>
#include <tr1/tuple>

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Enums/Dimensions.hpp"
#include "Exceptions/Exception.hpp"
#include "IoHdf5/Hdf5Types.hpp"
#include "IoHdf5/Hdf5File.hpp"

namespace GeoMHDiSCC {

namespace IoHdf5 {

   /**
    * @brief Interface to a general HDF5 file reader
    */
   class IHdf5Reader: public Hdf5File
   {
      public:
         /**
          * @brief Constructor
          *
          * @param name     Filename
          * @param ext      File extension
          * @param header   Header string of file
          * @param type     Type string of file
          * @param version  Version string of file
          */
         IHdf5Reader(std::string name, std::string ext, std::string header, std::string type, std::string version);

         /**
         * @brief Destructor
         */
         virtual ~IHdf5Reader();

         /**
          * @brief Initialise the file
          */
         void init();

         /**
          * @brief Read the content
          */

         virtual void read() = 0;

         /**
          * @brief Finalise the file
          */
         virtual void finalize();
         
      protected:
         /**
          * @brief Number of blocks to read
          */
         hsize_t  mBlock;

         /**
          * @brief Max collective IO read operations over all CPUs
          */
         int mCollIoRead;

         /**
          * @brief Set data parameters
          */
         virtual void setReadArguments() = 0;

         /**
          * @brief Open the file
          */
         void open();

         /**
          * @brief Close the file
          */
         void close();

         /**
          * @brief Check compatibility of opened file
          */
         void checkCompatibility();

         /**
          * @brief Read scalar dataset
          *
          * \warning In the MPI case only the IO node is going to read data
          *
          * @param loc     Base location in HDF5 file
          * @param dsname  HDF5 dataset name
          * @param val     Storage for scalar to read
          *
          * \tparam T type of the scalar
          */
         template <typename T> void readScalar(hid_t loc, const std::string dsname, T& val);

         /**
          * @brief Read Array dataset
          *
          * \warning In the MPI case only the IO node is going to read data
          *
          * @param loc     Base location in HDF5 file
          * @param dsname  HDF5 dataset name
          * @param arr     Storage for the data to read
          *
          * \tparam T type of the scalars
          */
         template <typename T> void readArray(hid_t loc, const std::string dsname, Eigen::Matrix<T, Eigen::Dynamic, 1>& arr);

         /**
          * @brief Read Matrix dataset
          *
          * \warning In the MPI case only the IO node is going to read data
          *
          * @param loc     Base location in HDF5 file
          * @param dsname  HDF5 dataset name
          * @param mat     Storage for the data to read
          *
          * \tparam T type of the scalars
          */
         template <typename T> void readMatrix(hid_t loc, const std::string dsname, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mat);

         /**
          * @brief Read an irregular field vector dataset
          *
          * @param loc     Base location in HDF5 file
          * @param dsname  HDF5 dataset name
          * @param storage Storage for the data to read (0 -> rows, 1 -> cols, 2 -> data pointer)
          *
          * \tparam T Type of the scalars
          */
         template <typename T> void readIrregularField(hid_t loc, const std::string dsname, std::vector<std::tr1::tuple<int, int, T *> >& storage);

         /**
          * @brief Read a regular field vector dataset
          *
          * @param loc     Base location in HDF5 file
          * @param dsname  HDF5 dataset name
          * @param storage Storage for the data to read (0 -> rows, 1 -> cols, 2 -> data pointer)
          *
          * \tparam T Type of the scalars
          */
         template <typename T> void readRegularField(hid_t loc, const std::string dsname, std::vector<std::tr1::tuple<int, int, T *> >& storage);

      private:
   };

   template <typename T> void IHdf5Reader::readScalar(hid_t loc, const std::string dsname, T& val)
   {
      // Open dataset
      hid_t dataset = H5Dopen(loc, dsname.c_str(), H5P_DEFAULT);

      // Read dataset with right data type
      H5Dread(dataset, Hdf5Types::type<T>(), H5S_ALL, H5S_ALL, H5P_DEFAULT, &val);

      // Close dataset
      H5Dclose(dataset);
   }

   template <typename T> void IHdf5Reader::readArray(hid_t loc, const std::string dsname, Eigen::Matrix<T, Eigen::Dynamic, 1>& arr)
   {
      // Set data type correctly
      hid_t type = Hdf5Types::type<T>();

      // Open dataset in file
      hid_t dataset = H5Dopen(loc, dsname.c_str(), H5P_DEFAULT);

      // Get file dataspace
      hid_t filespace = H5Dget_space(dataset);

      // Get rank of data set
      int rank = H5Sget_simple_extent_ndims(filespace);

      // Rank should be one
      if(rank != 1)
      {
         throw Exception("Rank of HDF5 dataset is not 1!");
      }

      // Read dimensions
      hsize_t dim;
      H5Sget_simple_extent_dims(filespace, &dim, NULL);

      // Supplied array should have right size (don't want to resize it)
      if(arr.size() != dim)
      {
         throw Exception("Provided storage for HDF5 data has wrong size!");
      }

      // memory dataspace 
      hid_t memspace = H5Screate_simple(1, &dim, NULL);

      // Read data from file
      H5Dread(dataset, type, memspace, filespace, H5P_DEFAULT, arr.data());
      
      // Close memspace
      H5Sclose(memspace);

      // Close filespace
      H5Sclose(filespace);

      // Close Dataset
      H5Dclose(dataset);
   }

   template <typename T> void IHdf5Reader::readMatrix(hid_t loc, const std::string dsname, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mat)
   {
      // Set data type correctly
      hid_t type = Hdf5Types::type<T>();

      // Open dataset in file
      hid_t dataset = H5Dopen(loc, dsname.c_str(), H5P_DEFAULT);

      // Get file dataspace
      hid_t filespace = H5Dget_space(dataset);

      // Get rank of data set
      int rank = H5Sget_simple_extent_ndims(filespace);

      // Rank should be one
      if(rank != 2)
      {
         throw Exception("Rank of HDF5 dataset is not 2!");
      }

      // Read dimensions
      hsize_t dims[2];
      H5Sget_simple_extent_dims(filespace, dims, NULL);

      // Supplied array should have right size (don't want to resize it)
      if(mat.cols() != dims[0] && mat.rows() != dims[1])
      {
         throw Exception("Provided storage for HDF5 data has wrong size!");
      }

      // memory dataspace 
      hid_t memspace = H5Screate_simple(2, dims, NULL);

      // Read data from file
      H5Dread(dataset, type, memspace, filespace, H5P_DEFAULT, mat.data());
      
      // Close memspace
      H5Sclose(memspace);

      // Close filespace
      H5Sclose(filespace);

      // Close Dataset
      H5Dclose(dataset);
   }

   template <typename T> void IHdf5Reader::readRegularField(hid_t loc, const std::string dsname, std::vector<std::tr1::tuple<int, int, T *> >& storage)
   {
      // Set data type correctly
      hid_t  type = Hdf5Types::type<T>();

      // Open dataset in file
      hid_t dataset = H5Dopen(loc, dsname.c_str(), H5P_DEFAULT);

      // Get file dataspace
      hid_t  filespace = H5Dget_space(dataset);

      // Get dimensionality of data set
      int nDims = H5Sget_simple_extent_ndims(filespace);

      // Read file dimensions
      hsize_t fDims[nDims];
      H5Sget_simple_extent_dims(filespace, fDims, NULL);

      // memory dataspace 
      hid_t  memspace;

      // Create offset storage
      hsize_t pOffset[nDims];
      // Fix offset for the fastest dimension (always full on each CPU)
      pOffset[nDims-1] = 0;

      // Compute size of the dataspace
      hsize_t dims[nDims];
      // Always read one block of the slowest dimension
      dims[0] = 1;
      // Fastest read dimension is minimum between file and memory sizes
      dims[nDims-1] = this->mBlock;

      // Compute size of the memory dataspace
      hsize_t iDims[2];

      // Storage for the memory offsets
      hsize_t memOffset[2];
      memOffset[0] = 0;
      memOffset[1] = 0;

      // Create dataset PList
      hid_t dsPList = this->datasetPList();

      //
      // First do the collective reads
      //
      for(int i = 0; i < this->mCollIoRead; ++i)
      {
         // Set full memory space
         iDims[0] = std::tr1::get<1>(storage.at(i));
         iDims[1] = std::tr1::get<0>(storage.at(i));
         memspace = H5Screate_simple(2, iDims, NULL);

         // Select memory hyperslabs (i.e. minimum between file and memory spaces)
         // Check that at list one index is stored in file
         if(this->mFileOffsets.at(i).size() > 0)
         {
            iDims[0] = std::min(std::tr1::get<1>(storage.at(i)), std::max(0,static_cast<int>(fDims[nDims-2]-this->mFileOffsets.at(i).at(1))));
         } else
         {
            iDims[0] = 0;
         }
         iDims[1] = this->mBlock;
         H5Sselect_hyperslab(memspace, H5S_SELECT_SET, memOffset, NULL, iDims, NULL);

         // Set offsets
         for(int j = 0; j < nDims -1; j++)
         {
            pOffset[j] = this->mFileOffsets.at(i).at(j);
         }

         // Select corresponding hyperslab
         dims[nDims-2] = iDims[0];
         H5Sselect_hyperslab(filespace, H5S_SELECT_SET, pOffset, NULL, dims, NULL);

         // Read file hyperslab into memory hyperslab
         H5Dread(dataset, type, memspace, filespace, dsPList, std::tr1::get<2>(storage.at(i)));

         // Reset hyperslab to whole dataset
         H5Sselect_all(filespace);

         // Close memory space
         H5Sclose(memspace);
      }

      //
      // Then do the independent reads
      //
      for(unsigned int i = this->mCollIoRead; i < this->mFileOffsets.size(); ++i)
      {
         // Set full memory space
         iDims[0] = std::tr1::get<1>(storage.at(i));
         iDims[1] = std::tr1::get<0>(storage.at(i));
         memspace = H5Screate_simple(2, iDims, NULL);

         // Select memory hyperslabs
         iDims[0] = std::min(std::tr1::get<1>(storage.at(i)), static_cast<int>(fDims[nDims-2]));
         iDims[1] = this->mBlock;
         H5Sselect_hyperslab(memspace, H5S_SELECT_SET, memOffset, NULL, iDims, NULL);

         // Set offsets
         for(int j = 0; j < nDims -1; j++)
         {
            pOffset[j] = this->mFileOffsets.at(i).at(j);
         }

         // Select corresponding hyperslab
         dims[nDims-2] = std::min(std::tr1::get<1>(storage.at(i)), static_cast<int>(fDims[nDims-2]));
         H5Sselect_hyperslab(filespace, H5S_SELECT_SET, pOffset, NULL, dims, NULL);

         // Read file hyperslab into memory hyperslab and update offset
         H5Dread(dataset, type, memspace, filespace, H5P_DEFAULT, std::tr1::get<2>(storage.at(i)));

         // Reset hyperslab to whole dataset
         H5Sselect_all(filespace);

         // Close memory space
         H5Sclose(memspace);

      }

      // Close file dataspace
      H5Sclose(filespace);

      // Close Dataset
      H5Dclose(dataset);

      this->freePList(dsPList);
   }

   template <typename T> void IHdf5Reader::readIrregularField(hid_t loc, const std::string dsname, std::vector<std::tr1::tuple<int, int, T *> >& storage)
   {
      // Set data type correctly
      hid_t  type = Hdf5Types::type<T>();

      // Open dataset in file
      hid_t dataset = H5Dopen(loc, dsname.c_str(), H5P_DEFAULT);

      // Get file dataspace
      hid_t  filespace = H5Dget_space(dataset);

      // memory dataspace 
      hid_t  memspace;

      // Storage for the memory offsets
      hsize_t memOffset[2];
      memOffset[0] = 0;
      memOffset[1] = 0;

      // Create offset storage
      hsize_t pOffset[2];
      pOffset[1] = 0;

      // Compute size of the memory dataspace
      hsize_t iDims[2];

      // Create dataset PList
      hid_t dsPList = this->datasetPList();

      //
      // First do the collective reads
      //
      for(int i = 0; i < this->mCollIoRead; ++i)
      {
         // Set full memory space
         iDims[0] = std::tr1::get<1>(storage.at(i));
         iDims[1] = std::tr1::get<0>(storage.at(i));
         memspace = H5Screate_simple(2, iDims, NULL);

         // Set size of read block (file resolution and memory resolution might be different)
         iDims[1] = this->mBlock;
         // Select memory hyperslabs
         H5Sselect_hyperslab(memspace, H5S_SELECT_SET, memOffset, NULL, iDims, NULL);

         // Create non regular hyperslab selection in file
         H5Sselect_none(filespace);
         iDims[0] = 1;
         for(unsigned int j =0; j < this->mFileOffsets.at(i).size(); ++j)
         {
            pOffset[0] = this->mFileOffsets.at(i).at(j);
            H5Sselect_hyperslab(filespace, H5S_SELECT_OR, pOffset, NULL, iDims, NULL);
         }

         // Read file hyperslab into memory hyperslab and update offset
         H5Dread(dataset, type, memspace, filespace, dsPList, std::tr1::get<2>(storage.at(i)));

         // Reset hyperslab to whole dataset
         H5Sselect_all(filespace);

         // Close memory space
         H5Sclose(memspace);
      }

      //
      // Then do the independent reads
      //
      for(unsigned int i = this->mCollIoRead; i < this->mFileOffsets.size(); ++i)
      {
         // Set full memory space
         iDims[0] = std::tr1::get<1>(storage.at(i));
         iDims[1] = std::tr1::get<0>(storage.at(i));
         memspace = H5Screate_simple(2, iDims, NULL);

         // Select memory hyperslabs
         iDims[1] = this->mBlock;
         H5Sselect_hyperslab(memspace, H5S_SELECT_SET, memOffset, NULL, iDims, NULL);

         // Create non regular hyperslab selection in file
         H5Sselect_none(filespace);
         iDims[0] = 1;
         for(unsigned int j =0; j < this->mFileOffsets.at(i).size(); ++j)
         {
            pOffset[0] = this->mFileOffsets.at(i).at(j);
            H5Sselect_hyperslab(filespace, H5S_SELECT_OR, pOffset, NULL, iDims, NULL);
         }

         // Read file hyperslab into memory hyperslab and update offset
         H5Dread(dataset, type, memspace, filespace, H5P_DEFAULT, std::tr1::get<2>(storage.at(i)));

         // Reset hyperslab to whole dataset
         H5Sselect_all(filespace);

         // Close memory space
         H5Sclose(memspace);
      }

      // Close file dataspace
      H5Sclose(filespace);

      // Close Dataset
      H5Dclose(dataset);

      this->freePList(dsPList);
   }

   /// Typedef for a shared pointer of a IHdf5Reader
   typedef SharedPtrMacro<IHdf5Reader> SharedIHdf5Reader;

}
}

#endif // IHDF5READER_HPP
