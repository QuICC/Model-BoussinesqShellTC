/** 
 * @file IAsciiNWriter.cpp
 * @brief Source of the interface to a numbered file ASCII writer
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 * @version 0.9.0
 * @date 2013-09-11
 */

// Configuration includes
//
#include "Framework/FrameworkMacro.h"

// System includes
//
#include <sstream>
#include <iomanip>

// External includes
//

// Class include
//
#include "IoAscii/IAsciiNWriter.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace IoAscii {

   const int IAsciiNWriter::msIDWidth = 4;

   IAsciiNWriter::IAsciiNWriter(std::string name, std::string ext, std::string header, std::string type, std::string version)
      : IAsciiWriter(name, ext, header, type, version), mCounter(0), mBaseName(name)
   {
   }

   IAsciiNWriter::~IAsciiNWriter()
   {
   }

   void IAsciiNWriter::preWrite()
   {
      // Check if the workflow allows IO to be performed
      if(FrameworkMacro::allowsIO())
      {
         // Update name with counter value
         this->updateName();

         // Create file
         this->open();

         // Add header information
         this->mFile << this->header() << std::endl;
         this->mFile << this->type() << std::endl;
         this->mFile << this->version() << std::endl;
      }
   }

   void IAsciiNWriter::postWrite()
   {
      // Check if the workflow allows IO to be performed
      if(FrameworkMacro::allowsIO())
      {
         // Close the file
         this->close();

         // Increment the file counter
         ++this->mCounter;
      }
   }

   void IAsciiNWriter::init()
   {
   }

   void IAsciiNWriter::finalize()
   {
   }

   void IAsciiNWriter::updateName()
   {
      // Create stringstream
      std::ostringstream   oss;

      // Create zerofilled string out of counter value
      oss << std::setfill('0') << std::setw(msIDWidth) << this->mCounter;

      // Reset filename including counter zerofilled value
      this->resetName(this->mBaseName + oss.str());
   }

}
}
