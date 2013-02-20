/** \file AsciiFile.cpp
 *  \brief Source of the general ASCII file implementation
 */

// System includes
//

// External includes
//

// Class include
//
#include "IoAscii/AsciiFile.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace IoAscii {

   AsciiFile::AsciiFile(std::string name, std::string ext, std::string header, std::string type, std::string version)
      : mName(name), mExt(ext), mHeader(header), mType(type), mVersion(version)
   {
   }

   std::string AsciiFile::filename()
   {
      return this->mName + this->mExt;
   }

   void AsciiFile::reset(std::string name)
   {
      this->mName = name;
   }

   std::string AsciiFile::name()
   {
      return this->mName;
   }

   std::string AsciiFile::extension()
   {
      return this->mExt;
   }

   std::string AsciiFile::header()
   {
      return AsciiFile::FILE_HEADER + this->mHeader;
   }

   std::string AsciiFile::type()
   {
      return AsciiFile::FILE_TYPE + this->mType;
   }

   std::string AsciiFile::version()
   {
      return AsciiFile::FILE_VERSION + this->mVersion;
   }

   const std::string AsciiFile::FILE_HEADER = "#";

   const std::string AsciiFile::FILE_TYPE = "#Type: ";

   const std::string AsciiFile::FILE_VERSION = "#Version: ";
}
}
