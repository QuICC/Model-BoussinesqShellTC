/** 
 * @file DirectBinaryWriter.cpp
 * @brief Source of the direct access binary writer
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
#include "IoBinary/DirectBinaryWriter.hpp"

// Project includes
//

namespace QuICC {

namespace IoBinary {

   DirectBinaryWriter::DirectBinaryWriter(std::string name, std::string ext, std::string header, std::string type, std::string version)
      : IBinaryWriter(name, ext, header, type, version)
   {
   }

   DirectBinaryWriter::~DirectBinaryWriter()
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

   void DirectBinaryWriter::finalize()
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
