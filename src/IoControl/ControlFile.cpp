/** 
 * @file ControlFile.cpp
 * @brief Source of the control file implementation
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// System includes
//

// External includes
//

// Class include
//
#include "IoControl/ControlFile.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace IoControl {

   ControlFile::ControlFile()
   {
   }

   ControlFile::~ControlFile()
   {
   }

   std::string ControlFile::filename() const
   {
      return ControlFile::FILE_NAME + ControlFile::FILE_EXTENSION;
   }

   std::string ControlFile::name() const
   {
      return ControlFile::FILE_NAME;
   }

   std::string ControlFile::extension() const
   {
      return ControlFile::FILE_EXTENSION;
   }

   std::string ControlFile::header() const
   {
      return ControlFile::FILE_HEADER;
   }

   std::string ControlFile::type() const
   {
      return ControlFile::FILE_TYPE;
   }

   std::string ControlFile::version() const
   {
      return ControlFile::FILE_VERSION;
   }

   const std::string ControlFile::FILE_NAME = "CONTROL";

   const std::string ControlFile::FILE_HEADER = "ControlInterface";

   const std::string ControlFile::FILE_TYPE = "ASCII";

   const std::string ControlFile::FILE_VERSION = "1.0";

   const std::string ControlFile::FILE_EXTENSION = "";
}
}
