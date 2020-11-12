/** 
 * @file IAsciiWriter.cpp
 * @brief Source of the interfact to a general ASCII writer
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
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
#include "IoAscii/IAsciiWriter.hpp"

// Project includes
//
#include "Exceptions/Exception.hpp"

namespace QuICC {

namespace IoAscii {

   const int IAsciiWriter::msIDWidth = 4;

   IAsciiWriter::IAsciiWriter(std::string name, std::string ext, std::string header, std::string type, std::string version, const IAsciiWriter::WriteMode mode)
      : AsciiFile(name, ext, header, type, version), mCounter(0), mBaseName(name), mMode(mode), mIsInitialized(false)
   {
   }

   IAsciiWriter::~IAsciiWriter()
   {
   }

   void IAsciiWriter::overwriteOutput()
   {
      if(this->mIsInitialized)
      {
         throw Exception("Cannot change output mode of " + this->filename() + " after initialization!");
      } else
      {
         this->mMode = IAsciiWriter::OVERWRITE;
      }
   }

   void IAsciiWriter::numberOutput()
   {
      if(this->mIsInitialized)
      {
         throw Exception("Cannot change output mode of " + this->filename() + " after initialization!");
      } else
      {
         this->mMode = IAsciiWriter::NUMBER;
      }
   }

   void IAsciiWriter::extendOutput()
   {
      if(this->mIsInitialized)
      {
         throw Exception("Cannot change output mode of " + this->filename() + " after initialization!");
      } else
      {
         this->mMode = IAsciiWriter::EXTEND;
      }
   }

   void IAsciiWriter::init()
   {
      this->mIsInitialized = true;

      if(this->mMode == IAsciiWriter::EXTEND)
      {
         this->create();
      }
   }

   void IAsciiWriter::finalize()
   {
      if(this->mMode == IAsciiWriter::EXTEND)
      {
         this->end();
      }
   }

   void IAsciiWriter::preWrite()
   {
      if(this->mMode == IAsciiWriter::OVERWRITE || this->mMode == IAsciiWriter::NUMBER)
      {
         this->create();
      }
   }

   void IAsciiWriter::postWrite()
   {
      if(this->mMode == IAsciiWriter::OVERWRITE || this->mMode == IAsciiWriter::NUMBER)
      {
         this->end();
      }
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

   void IAsciiWriter::openDebug()
   {
      // Get file handle
      this->mFile.open(this->filename().c_str());

      // Check that opening file was successful
      if(! this->mFile.is_open())
      {
         throw Exception("Couldn't open ASCII file " + this->filename() + "!");
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

   void IAsciiWriter::closeDebug()
   {
      this->mFile.close();
   }

   void IAsciiWriter::create()
   {
      // Check if the workflow allows IO to be performed
      if(FrameworkMacro::allowsIO())
      {
         if(this->mMode == IAsciiWriter::NUMBER)
         {
            this->updateName();
         }

         // Create file
         this->open();

         // Add header information
         this->mFile << this->header() << std::endl;
         this->mFile << this->type() << std::endl;
         this->mFile << this->version() << std::endl;
      }
   }

   void IAsciiWriter::end()
   {
      // Check if the workflow allows IO to be performed
      if(FrameworkMacro::allowsIO())
      {
         // Close the file
         this->close();

         if(this->mMode == IAsciiWriter::NUMBER)
         {
            // Increment the file counter
            ++this->mCounter;
         }
      }
   }

   void IAsciiWriter::updateName()
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
