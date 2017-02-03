/** 
 * @file IBinaryWriter.cpp
 * @brief Source of the general binary writer
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
#include "IoBinary/IBinaryWriter.hpp"

// Project includes
//
#include "Exceptions/Exception.hpp"

namespace QuICC {

namespace IoBinary {

   IBinaryWriter::IBinaryWriter(std::string name, std::string ext, std::string header, std::string type, std::string version)
      : BinaryFile(name, ext, header, type, version)
   {
   }

   IBinaryWriter::~IBinaryWriter()
   {
   }

   void IBinaryWriter::open()
   {
      // Check if the workflow allows IO to be performed
      if(FrameworkMacro::allowsIO())
      {
         // Get file handle
         this->mFile.open(this->filename().c_str(), std::ios::out | std::ios::binary);

         // Check that opening file was successful
         if(! this->mFile.is_open())
         {
            throw Exception("Couldn't open file " + this->filename() + "!");
         }
      }
   }

   void IBinaryWriter::close()
   {
      // Check if the workflow allows IO to be performed
      if(FrameworkMacro::allowsIO())
      {
         this->mFile.close();
      }
   }

}
}
