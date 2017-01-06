/** 
 * @file IAsciiEWriter.cpp
 * @brief Source of the interface to an "extending" ASCII writer
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
#include "IoAscii/IAsciiEWriter.hpp"

// Project includes
//

namespace QuICC {

namespace IoAscii {

   IAsciiEWriter::IAsciiEWriter(std::string name, std::string ext, std::string header, std::string type, std::string version)
      : IAsciiWriter(name, ext, header, type, version)
   {
   }

   IAsciiEWriter::~IAsciiEWriter()
   {
   }

   void IAsciiEWriter::init()
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

   void IAsciiEWriter::finalize()
   {
      // Check if the workflow allows IO to be performed
      if(FrameworkMacro::allowsIO())
      {
         // Close the file
         this->close();
      }
   }

   void IAsciiEWriter::preWrite()
   {
   }

   void IAsciiEWriter::postWrite()
   {
   }

}
}
