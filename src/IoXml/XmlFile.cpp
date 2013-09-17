/** 
 * @file XmlFile.cpp
 * @brief Source of the general XML file implementation
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
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
      : mXML(), mName(name), mExt(ext), mHeader(header), mType(type), mVersion(version)
   {
   }

   XmlFile::~XmlFile()
   {
   }

   std::string XmlFile::filename() const
   {
      return this->mName + this->mExt;
   }

   void XmlFile::resetName(std::string name)
   {
      this->mName = name;
   }

   std::string XmlFile::name() const
   {
      return this->mName;
   }

   std::string XmlFile::extension() const
   {
      return this->mExt;
   }

   std::string XmlFile::fileTag() const
   {
      return XmlFile::FILE_TAG;
   }

   std::string XmlFile::headerTag() const
   {
      return XmlFile::HEADER_TAG;
   }

   std::string XmlFile::header() const
   {
      return this->mHeader;
   }

   std::string XmlFile::typeTag() const
   {
      return XmlFile::TYPE_TAG;
   }

   std::string XmlFile::type() const
   {
      return this->mType;
   }

   std::string XmlFile::versionTag() const
   {
      return XmlFile::VERSION_TAG;
   }

   std::string XmlFile::version() const
   {
      return this->mVersion;
   }

   const std::string XmlFile::FILE_TAG = "file";

   const std::string XmlFile::HEADER_TAG = "header";

   const std::string XmlFile::TYPE_TAG = "type";

   const std::string XmlFile::VERSION_TAG = "version";
}
}
