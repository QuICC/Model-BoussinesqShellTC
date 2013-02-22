/** \file Hdf5File.hpp
 *  \brief Implementation of a general HDF5 file
 */

#ifndef HDF5FILE_HPP
#define HDF5FILE_HPP

// System includes
//
#include <vector>
#include <hdf5.h>

// External includes
//

// Project includes
//

namespace GeoMHDiSCC {

namespace IoHdf5 {

   /**
    * @brief Implementation of a general HDF5 file
    */
   class Hdf5File
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
         Hdf5File(std::string name, std::string ext, std::string header, std::string type, std::string version);

         /**
          * @brief Destructor
          */
         virtual ~Hdf5File() {};

         /**
          * @brief Get filename
          */
         std::string  filename() const;
         
      protected:
         /**
          * @brief Get the name
          */
         std::string  name() const;

         /**
          * @brief Get extension
          */
         std::string  extension() const;

         /**
          * @brief Get the header tag
          */
         std::string  headerTag() const;

         /**
          * @brief Get the header content
          */
         std::string  header() const;

         /**
          * @brief Get the type tag
          */
         std::string  typeTag() const;

         /**
          * @brief Get the type content
          */
         std::string  type() const;

         /**
          * @brief Get the version tag
          */
         std::string  versionTag() const;

         /**
          * @brief Get the version content
          */
         std::string  version() const;

         /**
          * @brief Reset name
          *
          * @param name New name
          */
         void resetName(std::string name);

         /**
          * @brief HDF5 file handle
          */
         hid_t  mFile;

         /**
          * @brief Size of the dimensions (HDF5 cols)
          */
         std::vector<hsize_t>  mFileDims;

         /**
          * @brief vector of offsets
          */
         std::vector< std::vector<hsize_t> > mFileOffsets;

         /**
          * @brief Get the dataset property list
          */
         hid_t datasetPList();

         /**
          * @brief Get the file property list
          */
         hid_t filePList();

         /**
          * @brief Free the dataset property list
          *
          * @param dsPList    Dataset property list
          */
         void freePList(hid_t dsPList);

      private:
         /**
          * @brief Header of file before header information
          */
         static const std::string HEADER_TAG;

         /**
          * @brief Header of file before optional information
          */
         static const std::string TYPE_TAG;

         /**
          * @brief Header of file before version information
          */
         static const std::string VERSION_TAG;

         /**
          * @brief Name of the file without extension
          */
         std::string mName;

         /**
          * @brief File extension
          */
         std::string mExt;

         /**
          * @brief Header of the file to check compatibility
          */
         std::string mHeader;

         /**
          * @brief Type of the file to check compatibility
          */
         std::string mType;

         /**
          * @brief Version of the file to check compatibility
          */
         std::string mVersion;
   };

}
}

#endif // HDF5FILE_HPP