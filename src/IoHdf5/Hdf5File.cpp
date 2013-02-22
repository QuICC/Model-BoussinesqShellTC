/** \file Hdf5File.cpp
 *  \brief Source of the implementation of a general HDF5 file
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "IoHdf5/Hdf5File.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace IoHdf5 {

   Hdf5File::Hdf5File(std::string name, std::string ext, std::string header, std::string type, std::string version)
      : mName(name), mExt(ext), mHeader(header), mType(type), mVersion(version), mFile()
   {
   }

   hid_t Hdf5File::filePList()
   {
      // 
      // Default options are fine for serial case
      //

      #ifdef GEOMHDISCC_MPI
         // Create file access property list
         hid_t fPList = H5Pcreate(H5P_FILE_ACCESS);

         // Create the MPI IO access property
         H5Pset_fapl_mpio(fPList, MPI_COMM_WORLD, MPI_INFO_NULL);
      #else
         hid_t fPList(H5P_DEFAULT);
      #endif // GEOMHDISCC_MPI

      return fPList;
   }

   hid_t Hdf5File::datasetPList()
   {
      // 
      // Default options are fine for serial case
      //

      #ifdef GEOMHDISCC_MPI
         // Create dataset transfer property list
         hid_t dsPList = H5Pcreate(H5P_DATASET_XFER);

         // Set the transfer property to collective IO
         H5Pset_dxpl_mpio(dsPList, H5FD_MPIO_COLLECTIVE);

         // ... or set the transfer property to independent IO
         //H5Pset_dxpl_mpio(dsPList, H5FD_MPIO_INDEPENDENT);
      #else
         hid_t dsPList(H5P_DEFAULT);
      #endif // GEOMHDISCC_MPI

      return dsPList;
   }

   void Hdf5File::freePList(hid_t dsPList)
   {
      // Let the HDF5 library close the dataset property list
      H5Pclose(dsPList);
   }

   std::string Hdf5File::filename()
   {
      return this->mName + this->mExt;
   }

   void Hdf5File::resetName(std::string name)
   {
      this->mName = name;
   }

   std::string Hdf5File::name()
   {
      return this->mName;
   }

   std::string Hdf5File::extension()
   {
      return this->mExt;
   }

   std::string Hdf5File::headerTag()
   {
      return Hdf5File::HEADER_TAG;
   }

   std::string Hdf5File::header()
   {
      return this->mHeader;
   }

   std::string Hdf5File::typeTag()
   {
      return Hdf5File::TYPE_TAG;
   }

   std::string Hdf5File::type()
   {
      return this->mType;
   }

   std::string Hdf5File::versionTag()
   {
      return Hdf5File::VERSION_TAG;
   }

   std::string Hdf5File::version()
   {
      return this->mVersion;
   }

   const std::string Hdf5File::HEADER_TAG = "header";

   const std::string Hdf5File::TYPE_TAG = "type";

   const std::string Hdf5File::VERSION_TAG = "version";
}
}