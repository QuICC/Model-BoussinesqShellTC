/** \file BinaryFile.cpp
 *  \brief Source of the general binary file implementation
 */

// System includes
//

// External includes
//

// Class include
//
#include "IoBinary/BinaryFile.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace IoBinary {

   BinaryFile::BinaryFile(std::string name, std::string ext, std::string header, std::string type, std::string version)
      : mName(name), mExt(ext), mHeader(header), mType(type), mVersion(version)
   {
   }

   std::string BinaryFile::filename()
   {
      return this->mName + this->mExt;
   }

   void BinaryFile::reset(std::string name)
   {
      this->mName = name;
   }

   std::string BinaryFile::name()
   {
      return this->mName;
   }

   std::string BinaryFile::extension()
   {
      return this->mExt;
   }

   std::string BinaryFile::header()
   {
      return BinaryFile::FILE_HEADER + this->mHeader;
   }

   std::string BinaryFile::type()
   {
      return BinaryFile::FILE_TYPE + this->mType;
   }

   std::string BinaryFile::version()
   {
      return BinaryFile::FILE_VERSION + this->mVersion;
   }

   const std::string BinaryFile::FILE_HEADER = "#";

   const std::string BinaryFile::FILE_TYPE = "#Type: ";

   const std::string BinaryFile::FILE_VERSION = "#Version: ";
}
}
