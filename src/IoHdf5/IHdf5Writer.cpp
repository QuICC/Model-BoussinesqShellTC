/** \file IHdf5Writer.cpp
 *  \brief Source of the general HDF5 writer
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "IoHdf5/IHdf5Writer.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace IoHdf5 {

   IHdf5Writer::IHdf5Writer(std::string name, std::string ext, std::string header, std::string type, std::string version)
      : Hdf5File(name, ext, header, type, version), mCollIOWrite(10000)
   {
   }

   void IHdf5Writer::open()
   {
      // Set Create Property list
      this->setFilePList(); 

      // Create file
      this->mFile = H5Fcreate(this->filename().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, this->filePList());
   }

   void IHdf5Writer::close()
   {
      // Free the file property list
      this->freeFilePList();

      // Close file
      H5Fclose(this->mFile);
   }

   void IHdf5Writer::createFileInfo()
   {
      // Select root of the file
      hid_t loc = H5Gopen(this->mFile, "/", H5P_DEFAULT);

      // Define some used variables
      hid_t attr;
      hid_t dspace;

      // Get the string data type
      hid_t type = H5Tcopy(H5T_C_S1);

      // Create header attribute
      dspace = H5Screate(H5S_SCALAR);
      H5Tset_size(type, this->header().size());
      attr = H5Acreate(loc, this->headerTag().c_str(), type, dspace, H5P_DEFAULT, H5P_DEFAULT);
      H5Awrite(attr, type, this->header().c_str());
      H5Aclose(attr);
      H5Sclose(dspace);

      // Create type attribute
      dspace = H5Screate(H5S_SCALAR);
      H5Tset_size (type, this->type().size());
      attr = H5Acreate(loc, this->typeTag().c_str(), type, dspace, H5P_DEFAULT, H5P_DEFAULT);
      H5Awrite(attr, type, this->type().c_str());
      H5Aclose(attr);
      H5Sclose(dspace);

      // Create version attribute
      dspace = H5Screate(H5S_SCALAR);
      H5Tset_size (type, this->version().size());
      attr = H5Acreate(loc, this->versionTag().c_str(), type, dspace, H5P_DEFAULT, H5P_DEFAULT);
      H5Awrite(attr, type, this->version().c_str());
      H5Aclose(attr);
      H5Sclose(dspace);

      H5Gclose(loc);
   }
}
}