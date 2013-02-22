/** \file IAsciiWriter.cpp
 *  \brief Source of the interfact to a general ASCII writer
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
#include "IoAscii/IAsciiWriter.hpp"

// Project includes
//
#include "Exception/Exception.hpp"

namespace GeoMHDiSCC {

namespace IoAscii {

   IAsciiWriter::IAsciiWriter(std::string name, std::string ext, std::string header, std::string type, std::string version)
      : AsciiFile(name, ext, header, type, version)
   {
   }

   void IAsciiWriter::open()
   {
      // Check if the workflow allows IO to be performed
      if(FrameworkMacro::allowsIO())
      {
         // Get file handle
         this->mFile.open(this->filename().c_str());

         // Check that opening file was successful
         if(! this->mFile.is_open())
         {
            throw Exception("Couldn't open ASCII file " + this->filename() + "!");
         }
      }
   }

   void IAsciiWriter::close()
   {
      // Check if the workflow allows IO to be performed
      if(FrameworkMacro::allowsIO())
      {
         this->mFile.close();
      }
   }

}
}