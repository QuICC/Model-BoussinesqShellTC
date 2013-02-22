/** \file IAsciiRWriter.cpp
 *  \brief Source of the interface to an overwriting ASCII writer
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
#include "IoAscii/IAsciiRWriter.hpp"

// Project includes
//

namespace GeoMDHiSCC {

namespace IoAscii {

   IAsciiRWriter::IAsciiRWriter(std::string name, std::string ext, std::string header, std::string type, std::string version)
      : IAsciiWriter(name, ext, header, type, version)
   {
   }

   void IAsciiRWriter::init()
   {
   }

   void IAsciiRWriter::finalise()
   {
   }

   void IAsciiRWriter::preWrite()
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

   void IAsciiRWriter::postWrite()
   {
      // Check if the workflow allows IO to be performed
      if(FrameworkMacro::allowsIO())
      {
         // Close the file
         this->close();
      }
   }

}
}