/** \file XmlFile.cpp
 *  \brief Source of the general XML file implementation
 */

// System includes
//

// External includes
//

// Class include
//
#include "IoXml/XmlFile.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace IoXml {

   XmlFile::XmlFile(std::string name, std::string ext, std::string header, std::string type, std::string version)
      : mName(name), mExt(ext), mHeader(header), mType(type), mVersion(version), mXML()
   {
   }

   std::string XmlFile::filename()
   {
      return this->mName + this->mExt;
   }

   void XmlFile::reset(std::string name)
   {
      this->mName = name;
   }

   std::string XmlFile::name()
   {
      return this->mName;
   }

   std::string XmlFile::extension()
   {
      return this->mExt;
   }

   std::string XmlFile::fileTag()
   {
      return XmlFile::FILE_TAG;
   }

   std::string XmlFile::headerTag()
   {
      return XmlFile::HEADER_TAG;
   }

   std::string XmlFile::header()
   {
      return this->mHeader;
   }

   std::string XmlFile::typeTag()
   {
      return XmlFile::TYPE_TAG;
   }

   std::string XmlFile::type()
   {
      return this->mType;
   }

   std::string XmlFile::versionTag()
   {
      return XmlFile::VERSION_TAG;
   }

   std::string XmlFile::version()
   {
      return this->mVersion;
   }

   const std::string XmlFile::FILE_TAG = "file";

   const std::string XmlFile::HEADER_TAG = "header";

   const std::string XmlFile::TYPE_TAG = "type";

   const std::string XmlFile::VERSION_TAG = "version";
}
}
