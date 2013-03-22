/** \file XmlFile.hpp
 *  \brief Implementation of a general XML file
 */

#ifndef XMLFILE_HPP
#define XMLFILE_HPP

// System includes
//
#include <string>

// External includes
//
#include <rapidxml.hpp>

// Project includes
//

namespace GeoMHDiSCC {

namespace IoXml {

   /**
    * @brief Implementation of a general XML file
    */
   class XmlFile
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
         XmlFile(std::string name, std::string ext, std::string header, std::string type, std::string version);

         /**
         * @brief Destructor
         */
         virtual ~XmlFile();

         /**
          * @brief Get filename
          */
         std::string  filename() const;
         
      protected:
         /**
          * @brief XML interface
          */
         rapidxml::xml_document<>   mXML;

         /**
          * @brief Get the name
          */
         std::string  name() const;

         /**
          * @brief Get extension
          */
         std::string  extension() const;

         /**
          * @brief Get file meta tag
          */
         std::string fileTag() const;

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

      private:
         /**
          * @brief XML tag for meta file information
          */
         static const std::string FILE_TAG;

         /**
          * @brief XML tag for file header information
          */
         static const std::string HEADER_TAG;

         /**
          * @brief XML tag for file type information
          */
         static const std::string TYPE_TAG;

         /**
          * @brief XML tag for file version information
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

#endif // XMLFILE_HPP
