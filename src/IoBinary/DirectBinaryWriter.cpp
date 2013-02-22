/** \file DirectBinaryWriter.cpp
 *  \brief Source of the direct access binary writer
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
#include "IoBinary/DirectBinaryWriter.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace IoBinary {

   DirectBinaryWriter::DirectBinaryWriter(std::string name, std::string ext, std::string header, std::string type, std::string version)
      : IBinaryWriter(name, ext, header, type, version)
   {
   }

   void DirectBinaryWriter::init()
   {
      // Check if the workflow allows IO to be performed
      if(FrameworkMacro::allowsIO())
      {
         // Create file
         this->open();
      }
   }

   void DirectBinaryWriter::finalise()
   {
      // Check if the workflow allows IO to be performed
      if(FrameworkMacro::allowsIO())
      {
         // Close the file
         this->close();
      }
   }

   std::ofstream& DirectBinaryWriter::file()
   {
      return this->mFile;
   }

}
}