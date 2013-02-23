/** \file IHdf5Writer.hpp
 *  \brief Interface to a general HDF5 writer
 */

#ifndef IHDF5WRITER_HPP
#define IHDF5WRITER_HPP

// Configuration includes
//
#include "Framework/FrameworkMacro.h"
#include "SmartPointers/SharedPtrMacro.h"

// System includes
//
#include <string>

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "IoHdf5/Hdf5Types.hpp"
#include "IoHdf5/Hdf5File.hpp"

namespace GeoMHDiSCC {

namespace IoHdf5 {

   /**
    * @brief Interface to a general HDF5 file writer
    */
   class IHdf5Writer: public Hdf5File
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
         IHdf5Writer(std::string name, std::string ext, std::string header, std::string type, std::string version);

         /**
          * @brief Destructor
          */
         virtual ~IHdf5Writer();

         /**
          * @brief Initialise the file
          */
         virtual void init() = 0;

         /**
          * @brief Write the content
          */
         virtual void write() = 0;

         /**
          * @brief Finalise the file
          */
         virtual void finalize() = 0;
         
      protected:
         /**
          * @brief Max collective IO Write operations over all CPUs
          */
         int mCollIOWrite;

         /**
          * @brief Open the file
          */
         void open();

         /**
          * @brief Close the file
          */
         void close();

         /**
          * @brief Create the file info (header, type and version)
          */
         void createFileInfo();

         /**
          * @brief Write scalar dataset
          *
          * \warning In the MPI case only the IO node is going to write data
          *
          * @param loc     Base location in HDF5 file
          * @param dsname  HDF5 dataset name
          * @param val     Scalar to write
          *
          * \tparam T Type of the scalar
          */
         template <typename T> void writeScalar(hid_t loc, const std::string dsname, const T val);

         /**
          * @brief Write Array dataset
          *
          * \warning In the MPI case only the IO node is going to write data
          *
          * @param loc     Base location in HDF5 file
          * @param dsname  HDF5 dataset name
          * @param arr     Array of values to write
          *
          * \tparam T Type of the scalars
          */
         template <typename T> void writeArray(hid_t loc, const std::string dsname, const Eigen::Matrix<T, Eigen::Dynamic, 1>& arr);

         /**
          * @brief Write Matrix dataset
          *
          * \warning In the MPI case only the IO node is going to write data
          *
          * @param loc     Base location in HDF5 file
          * @param dsname  HDF5 dataset name
          * @param mat     Matrix of values to write
          *
          * \tparam T Type of the scalars
          */
         template <typename T> void writeMatrix(hid_t loc, const std::string dsname, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mat);

         /**
          * @brief Write an irregular field vector dataset
          *
          * For irregular field data the values gets stored in a 2D array
          *
          * \warning This method does not work for 1D data
          *
          * @param loc     Base location in HDF5 file
          * @param dsname  HDF5 dataset name
          * @param storage Special data to write
          *
          * \tparam T Type of the scalars
          */
         template <typename T> void writeIrregularField(hid_t loc, const std::string dsname, const std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >& storage);

         /**
          * @brief Write an irregular field vector dataset
          *
          * For irregular field data the values gets stored in a 2D array
          *
          * \warning This method does not work for 1D data
          *
          * @param loc     Base location in HDF5 file
          * @param dsname  HDF5 dataset name
          * @param storage Special data to write
          *
          * \tparam T Type of the scalars
          */
         template <typename T> void writeIrregularField(hid_t loc, const std::string dsname, const std::vector<Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > >& storage);

         /**
          * @brief Write a regular field vector dataset
          *
          * For regular field data the values gets stored in a array of the same
          * dimensionality as the data i.e. 3D in 3D array, 2D data in 2D array
          *
          * \warning This method does not work for 1D data
          *
          * @param loc     Base location in HDF5 file
          * @param dsname  HDF5 dataset name
          * @param storage Special data to write
          *
          * \tparam T Type of the scalars
          */
         template <typename T> void writeRegularField(hid_t loc, const std::string dsname, const std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >& storage);

         /**
          * @brief Write a regular field vector dataset (on sliced map)
          *
          * For regular field data the values gets stored in a array of the same
          * dimensionality as the data i.e. 3D in 3D array, 2D data in 2D array
          *
          * \warning This method does not work for 1D data
          *
          * @param loc     Base location in HDF5 file
          * @param dsname  HDF5 dataset name
          * @param storage Special data to write
          *
          * \tparam T Type of the scalars
          */
         template <typename T> void writeRegularField(hid_t loc, const std::string dsname, const std::vector<Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > >& storage);

      private:
   };

   template <typename T> void IHdf5Writer::writeScalar(hid_t loc, const std::string dsname, const T val)
   {
      // Get correct data type
      hid_t type = Hdf5Types::type<T>();

      // Create scalar dataspace
      hid_t dspace = H5Screate(H5S_SCALAR);

      // Create data set
      hid_t dataset = H5Dcreate(loc, dsname.c_str(), type, dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      // Write scalar value
      if(FrameworkMacro::allowsIO())
      {
         H5Dwrite(dataset, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, &val);
      }

      // close dataspace
      H5Sclose(dspace);

      // close dataset
      H5Dclose(dataset);
   }

   template <typename T> void IHdf5Writer::writeArray(hid_t loc, const std::string dsname, const Eigen::Matrix<T, Eigen::Dynamic, 1>& arr)
   {
      // Set data type correctly
      hid_t type = Hdf5Types::type<T>();

      // Compute size of the dataspace
      hsize_t dims = arr.size();

      // Create file dataspace
      hid_t filespace = H5Screate_simple(1, &dims, NULL);

      // Create dataset in file
      hid_t dataset = H5Dcreate(loc, dsname.c_str(), type, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      // Set memory space
      hid_t memspace = H5Screate_simple(1, &dims, NULL);

      // Write memory
      if(FrameworkMacro::allowsIO())
      {
         H5Dwrite(dataset, type, memspace, filespace, H5P_DEFAULT, arr.data());
      }
      
      // Close memspace
      H5Sclose(memspace);

      // Close filespace
      H5Sclose(filespace);

      // Close Dataset
      H5Dclose(dataset);
   }

   template <typename T> void IHdf5Writer::writeMatrix(hid_t loc, const std::string dsname, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mat)
   {
      // Set data type correctly
      hid_t type = Hdf5Types::type<T>();

      // Compute size of the dataspace
      hsize_t dims[2];
      dims[0] = mat.cols();
      dims[1] = mat.rows();

      // Create file dataspace
      hid_t filespace = H5Screate_simple(2, dims, NULL);

      // Create dataset in file
      hid_t dataset = H5Dcreate(loc, dsname.c_str(), type, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      // Set memory space
      hid_t memspace = H5Screate_simple(2, dims, NULL);

      // Write memory
      if(FrameworkMacro::allowsIO())
      {
         H5Dwrite(dataset, type, memspace, filespace, H5P_DEFAULT, mat.data());
      }
      
      // Close memspace
      H5Sclose(memspace);

      // Close filespace
      H5Sclose(filespace);

      // Close Dataset
      H5Dclose(dataset);
   }

   template <typename T> void IHdf5Writer::writeRegularField(hid_t loc, const std::string dsname, const std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >& storage)
   {
      // Get dimensionality of file data
      int nDims = this->mFileDims.size();

      // Set data type correctly
      hid_t type = Hdf5Types::type<T>();

      // Create file dataspace
      hid_t filespace = H5Screate_simple(nDims, &this->mFileDims.front(), NULL);

      // Create dataset in file
      hid_t dataset = H5Dcreate(loc, dsname.c_str(), type, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      // Memory dataspace 
      hid_t  memspace;

      // Create offset storage
      hsize_t pOffsets[nDims];
      pOffsets[nDims-1] = 0;

      // Compute size of the dataspace
      hsize_t dims[nDims];
      dims[0] = 1;
      dims[nDims-1] = this->mFileDims.back();

      // Compute size of the memory dataspace
      hsize_t iDims[2];
      iDims[1] = this->mFileDims.back();

      // Set PList to parallel access
      hid_t dsPList = this->datasetPList();

      // Loop over matrices to store in file
      for(unsigned int i = 0; i < storage.size() ; ++i)
      {
         // Set memory space
         iDims[0] = storage.at(i).cols();
         memspace = H5Screate_simple(2, iDims, NULL);

         // Set offsets
         for(int j = 0; j < nDims -1; j++)
         {
            pOffsets[j] = this->mFileOffsets.at(i).at(j);
         }

         // Select corresponding hyperslab
         dims[nDims-2] = storage.at(i).cols();
         H5Sselect_hyperslab(filespace, H5S_SELECT_SET, pOffsets, NULL, dims, NULL);

         // Write memory into hyperslab
         H5Dwrite(dataset, type, memspace, filespace, dsPList, storage.at(i).data());

         // Reset hyperslab to whole dataset
         H5Sselect_all(filespace);

         // Close memory space
         H5Sclose(memspace);
      }

      //
      // Add some ZERO IO calls to allow for collective writes
      //

      // Set zero memory space
      iDims[0] = 1;
      memspace = H5Screate_simple(1, iDims, NULL);
      H5Sselect_none(memspace);

      // Select zero file space
      H5Sselect_none(filespace);

      // Zero write IO
      for(unsigned int i = 0; i < (this->mCollIOWrite - storage.size()) ; ++i)
      {
         // Write memory into hyperslab
         H5Dwrite(dataset, type, memspace, filespace, dsPList, storage.at(0).data());
      }

      // Close memspace
      H5Sclose(memspace);

      // Close Dataset
      H5Dclose(dataset);

      // Close filespace
      H5Sclose(filespace);

      // Free the dataset PList
      this->freePList(dsPList);
   }

   template <typename T> void IHdf5Writer::writeRegularField(hid_t loc, const std::string dsname, const std::vector<Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > >& storageMap)
   {
      std::vector<Matrix> storage;
      for(int i = 0; i < storageMap.size(); i++)
      {
         storage.push_back(Matrix(storageMap.at(i).rows(), storageMap.at(i).rows()));
         storage.back() = storageMap.at(i);
      }

      // Get dimensionality of file data
      int nDims = this->mFileDims.size();

      // Set data type correctly
      hid_t type = Hdf5Types::type<T>();

      // Create file dataspace
      hid_t filespace = H5Screate_simple(nDims, &this->mFileDims.front(), NULL);

      // Create dataset in file
      hid_t dataset = H5Dcreate(loc, dsname.c_str(), type, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      // Memory dataspace 
      hid_t  memspace;

      // Create offset storage
      hsize_t pOffsets[nDims];
      pOffsets[nDims-1] = 0;

      // Compute size of the dataspace
      hsize_t dims[nDims];
      dims[0] = 1;
      dims[nDims-1] = this->mFileDims.back();

      // Compute size of the memory dataspace
      hsize_t iDims[2];
      iDims[1] = this->mFileDims.back();

      // Set PList to parallel access
      hid_t dsPList = this->datasetPList();

      // Loop over matrices to store in file
      for(unsigned int i = 0; i < storage.size() ; ++i)
      {
         // Set memory space
         iDims[0] = storage.at(i).cols();
         memspace = H5Screate_simple(2, iDims, NULL);

         // Set offsets
         for(int j = 0; j < nDims -1; j++)
         {
            pOffsets[j] = this->mFileOffsets.at(i).at(j);
         }

         // Select corresponding hyperslab
         dims[nDims-2] = storage.at(i).cols();
         H5Sselect_hyperslab(filespace, H5S_SELECT_SET, pOffsets, NULL, dims, NULL);

         // Write memory into hyperslab
         H5Dwrite(dataset, type, memspace, filespace, dsPList, storage.at(i).data());

         // Reset hyperslab to whole dataset
         H5Sselect_all(filespace);

         // Close memory space
         H5Sclose(memspace);
      }

      //
      // Add some ZERO IO calls to allow for collective writes
      //

      // Set zero memory space
      iDims[0] = 1;
      memspace = H5Screate_simple(1, iDims, NULL);
      H5Sselect_none(memspace);

      // Select zero file space
      H5Sselect_none(filespace);

      // Zero write IO
      for(unsigned int i = 0; i < (this->mCollIOWrite - storage.size()) ; ++i)
      {
         // Write memory into hyperslab
         H5Dwrite(dataset, type, memspace, filespace, dsPList, storage.at(0).data());
      }

      // Close memspace
      H5Sclose(memspace);

      // Close Dataset
      H5Dclose(dataset);

      // Close filespace
      H5Sclose(filespace);

      // Free the dataset PList
      this->freePList(dsPList);
   }

   template <typename T> void IHdf5Writer::writeIrregularField(hid_t loc, const std::string dsname, const std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >& storage)
   {
      // Get dimensionality of file data
      int nDims = this->mFileDims.size();

      // Set data type correctly
      hid_t type = Hdf5Types::type<T>();

      // Create file dataspace
      hid_t filespace = H5Screate_simple(nDims, &(this->mFileDims.front()), NULL);

      // Create dataset in file
      hid_t dataset = H5Dcreate(loc, dsname.c_str(), type, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      // Memory dataspace 
      hid_t  memspace;

      // Create offset storage
      hsize_t pOffset[2];
      pOffset[1] = 0;

      // Compute size of the memory dataspace
      hsize_t iDims[2];
      iDims[1] = this->mFileDims.back();

      // Set PList to parallel access
      hid_t dsPList = this->datasetPList();

      // Loop over matrices to store in file
      for(unsigned int i = 0; i < storage.size() ; ++i)
      {
         // Set memory space
         iDims[0] = storage.at(i).cols();
         memspace = H5Screate_simple(2, iDims, NULL);

         // Create non regular hyperslab selection in file
         H5Sselect_none(filespace);
         iDims[0] = 1;
         for(unsigned int j =0; j < this->mFileOffsets.at(i).size(); ++j)
         {
            pOffset[0] = this->mFileOffsets.at(i).at(j);
            H5Sselect_hyperslab(filespace, H5S_SELECT_OR, pOffset, NULL, iDims, NULL);
         }

         // Write memory into hyperslab
         H5Dwrite(dataset, type, memspace, filespace, dsPList, storage.at(i).data());

         // Reset hyperslab to whole dataset
         H5Sselect_all(filespace);

         // Close memory space
         H5Sclose(memspace);
      }

      //
      // Add some ZERO IO calls to allow for collective writes
      //
      
      // Create zero memory space
      iDims[0] = 1;
      memspace = H5Screate_simple(1, iDims, NULL);
      H5Sselect_none(memspace);

      // Create zero file selection space
      H5Sselect_none(filespace);

      // Zero write IO
      for(unsigned int i = 0; i < (this->mCollIOWrite - storage.size()); ++i)
      {
         // Write memory into hyperslab
         H5Dwrite(dataset, type, memspace, filespace, dsPList, storage.at(0).data());
      }

      // close memspace
      H5Sclose(memspace);

      // Close Dataset
      H5Dclose(dataset);

      // Close filespace
      H5Sclose(filespace);

      // Free the dataset PList
      this->freePList(dsPList);
   }

   template <typename T> void IHdf5Writer::writeIrregularField(hid_t loc, const std::string dsname, const std::vector<Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > >& storage)
   {
      // Get dimensionality of file data
      int nDims = this->mFileDims.size();

      // Set data type correctly
      hid_t type = Hdf5Types::type<T>();

      // Create file dataspace
      hid_t filespace = H5Screate_simple(nDims, &(this->mFileDims.front()), NULL);

      // Create dataset in file
      hid_t dataset = H5Dcreate(loc, dsname.c_str(), type, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      // Memory dataspace 
      hid_t  memspace;

      // Create offset storage
      hsize_t pOffset[2];
      pOffset[1] = 0;

      // Compute size of the memory dataspace
      hsize_t iDims[2];
      iDims[1] = this->mFileDims.back();

      // Set PList to parallel access
      hid_t dsPList = this->datasetPList();

      // Loop over matrices to store in file
      for(unsigned int i = 0; i < storage.size() ; ++i)
      {
         // Set memory space
         iDims[0] = storage.at(i).cols();
         memspace = H5Screate_simple(2, iDims, NULL);

         // Create non regular hyperslab selection in file
         H5Sselect_none(filespace);
         iDims[0] = 1;
         for(unsigned int j =0; j < this->mFileOffsets.at(i).size(); ++j)
         {
            pOffset[0] = this->mFileOffsets.at(i).at(j);
            H5Sselect_hyperslab(filespace, H5S_SELECT_OR, pOffset, NULL, iDims, NULL);
         }

         // Write memory into hyperslab
         H5Dwrite(dataset, type, memspace, filespace, dsPList, storage.at(i).data());

         // Reset hyperslab to whole dataset
         H5Sselect_all(filespace);

         // Close memory space
         H5Sclose(memspace);
      }

      //
      // Add some ZERO IO calls to allow for collective writes
      //
      
      // Create zero memory space
      iDims[0] = 1;
      memspace = H5Screate_simple(1, iDims, NULL);
      H5Sselect_none(memspace);

      // Create zero file selection space
      H5Sselect_none(filespace);

      // Zero write IO
      for(unsigned int i = 0; i < (this->mCollIOWrite - storage.size()); ++i)
      {
         // Write memory into hyperslab
         H5Dwrite(dataset, type, memspace, filespace, dsPList, storage.at(0).data());
      }

      // close memspace
      H5Sclose(memspace);

      // Close Dataset
      H5Dclose(dataset);

      // Close filespace
      H5Sclose(filespace);

      // Free the dataset PList
      this->freePList(dsPList);
   }

   /// Typedef for a shared pointer of a IHdf5Writer
   typedef SharedPtrMacro<IHdf5Writer> SharedIHdf5Writer;

}
}

#endif // IHDF5WRITER_HPP
