/** 
 * @file DirectAsciiWriter.cpp
 * @brief Source of the direct access ASCII writer
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Configuration includes
//
#include "Framework/FrameworkMacro.h"

// System includes
//

// External includes
//

// Class include
//
#include "IoAscii/DirectAsciiWriter.hpp"

// Project includes
//

namespace QuICC {

namespace IoAscii {

   DirectAsciiWriter::DirectAsciiWriter(std::string name, std::string ext, std::string header, std::string type, std::string version)
      : IAsciiWriter(name, ext, header, type, version)
   {
   }

   DirectAsciiWriter::~DirectAsciiWriter()
   {
   }

   void DirectAsciiWriter::init()
   {
      // Check if the workflow allows IO to be performed
      if(FrameworkMacro::allowsIO())
      {
         // Create file
         this->open();

         // Add header information
         this->mFile << this->header() << std::endl;
         this->mFile << this->type() << std::endl;
         this->mFile << this->version() << std::endl;
      }
   }

   void DirectAsciiWriter::initDebug()
   {
      // Create file
      this->openDebug();

      // Add header information
      this->mFile << this->header() << std::endl;
      this->mFile << this->type() << std::endl;
      this->mFile << this->version() << std::endl;
   }

   void DirectAsciiWriter::initNoHeader()
   {
      // Check if the workflow allows IO to be performed
      if(FrameworkMacro::allowsIO())
      {
         // Create file
         this->open();
      }
   }

   void DirectAsciiWriter::finalize()
   {
      // Check if the workflow allows IO to be performed
      if(FrameworkMacro::allowsIO())
      {
         // Close the file
         this->close();
      }
   }

   void DirectAsciiWriter::finalizeDebug()
   {
      // Close the file
      this->closeDebug();
   }

   std::ofstream& DirectAsciiWriter::file()
   {
      return this->mFile;
   }

}
}
