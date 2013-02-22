/** \file IHdf5Reader.cpp
 *  \brief Source of the HDF5 reader implementation
 */

// System includes
//
#include <sstream>

// External includes
//

// Class include
//
#include "IoHdf5/IHdf5Reader.hpp"

// Project includes
//
#include "Exceptions/Exception.hpp"

namespace GeoMHDiSCC {

namespace IoHdf5 {

   IHdf5Reader::IHdf5Reader(std::string name, std::string ext, std::string header, std::string type, std::string version)
      : Hdf5File(name, ext, header, type, version), mBlock(0), mCollIORead(0)
   {
   }

   void IHdf5Reader::init()
   {
      // Open the file
      this->open();

      // Check file compatibility
      this->checkCompatibility();
   }

   void IHdf5Reader::open()
   {
      // Open file
      this->mFile = H5Fopen(this->filename().c_str(), H5F_ACC_RDONLY, this->filePList());

      // Check for successfully opened file
      if(this->mFile < 0)
      {
         throw Exception("Failed to open HDF5 file " + this->filename() + "!");
      }
   }

   void IHdf5Reader::finalise()
   {
      this->close();
   }

   void IHdf5Reader::close()
   {
      hid_t fPList = H5Fget_access_plist(this->mFile);

      // Close file
      H5Fclose(this->mFile);

      // Free the property list
      this->freePList(fPList);
   }

   void IHdf5Reader::checkCompatibility()
   {
      // Open the root of the file
      hid_t loc = H5Gopen(this->mFile, "/", H5P_DEFAULT);

      // Some useful variables
      hid_t type, ftype;
      hid_t  attr;
      char *cStr;
      size_t size;

      //
      // Read header attribute
      //

      // Open the header attribute
      attr = H5Aopen(loc, this->headerTag().c_str(), H5P_DEFAULT);
      // Get type
      ftype = H5Aget_type(attr);
      // Get string size
      size = H5Tget_size(ftype);
      // Create string type
      type = H5Tcopy (H5T_C_S1);
      // Set string length
      H5Tset_size(type, size);
      // Alloate memory
      cStr = new char[size];
      // Read attribute
      H5Aread(attr, type, cStr);
      // Create string
      std::string fileHeader(cStr);
      // Check compatibility
      if(this->header() != fileHeader)
      {
         throw Exception("Wrong HDF5 file header!");
      }
      delete[] cStr;
      H5Aclose(attr);

      //
      // Read type attribute
      //

      // Open the type attribute
      attr = H5Aopen(loc, this->typeTag().c_str(), H5P_DEFAULT);
      // Get type
      ftype = H5Aget_type(attr);
      // Get string size
      size = H5Tget_size(ftype);
      // Create string type
      type = H5Tcopy (H5T_C_S1);
      // Set string length
      H5Tset_size(type, size);
      // Alloate memory
      cStr = new char[size];
      // Read attribute
      H5Aread(attr, type, cStr);
      // Create string
      std::string fileType(cStr);
      // Check compatibility
      if(this->type() != fileType)
      {
         throw Exception("Wrong HDF5 file type!");
      }
      delete[] cStr;
      H5Aclose(attr);

      //
      // Read version attribute
      //

      // Open the version attribute
      attr = H5Aopen(loc, this->versionTag().c_str(), H5P_DEFAULT);
      // Get type
      ftype = H5Aget_type(attr);
      // Get string size
      size = H5Tget_size(ftype);
      // Create string type
      type = H5Tcopy (H5T_C_S1);
      // Set string length
      H5Tset_size(type, size);
      // Alloate memory
      cStr = new char[size];
      // Read attribute
      H5Aread(attr, type, cStr);
      // Create string
      std::string fileVersion(cStr);
      // Check compatibility
      if(this->version() != fileVersion)
      {
         throw Exception("Wrong HDF5 file version!");
      }
      delete[] cStr;
      H5Aclose(attr);

      H5Gclose(loc);
   }
}
}
