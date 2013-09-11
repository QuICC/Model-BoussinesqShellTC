/** 
 * @file IAsciiReader.cpp
 * @brief Source of the interface to a general ASCII reader
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 * @version 0.9.0
 * @date 2013-09-11
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
#include "IoAscii/IAsciiReader.hpp"

// Project includes
//
#include "Exceptions/Exception.hpp"

namespace GeoMHDiSCC {

namespace IoAscii {

   IAsciiReader::IAsciiReader(std::string name, std::string ext, std::string header, std::string type, std::string version)
      : AsciiFile(name, ext, header, type, version)
   {
   }

   IAsciiReader::~IAsciiReader()
   {
   }

   void IAsciiReader::init()
   {
      // Check if the workflow allows IO to be performed
      if(FrameworkMacro::allowsIO())
      {
         // open the file
         this->open();

         // Check that the file is compatible with the reader
         this->checkCompatibility();
      }
   }

   void IAsciiReader::open()
   {
      // Check if the workflow allows IO to be performed
      if(FrameworkMacro::allowsIO())
      {
         // Get handle to file
         this->mFile.open(this->filename().c_str());

         // Check that opening was a success
         if(! this->mFile.is_open())
         {
            throw Exception("Couldn't open ASCII file " + this->filename() + "!");
         }
      }
   }

   void IAsciiReader::checkCompatibility()
   {
      // Check if the workflow allows IO to be performed
      if(FrameworkMacro::allowsIO())
      {
         std::string fileHead;
         std::string fileType;
         std::string fileVers;

         // Read header if present
         getline(this->mFile, fileHead);

         // Check reading header was successful
         if(this->mFile.good())
         {
            // Check header
            if(fileHead.compare(this->header()) == 0)
            {
               // Read type if present
               getline(this->mFile, fileType);

               // Check reading type was successful
               if(this->mFile.good())
               {
                  // Check type
                  if(fileType.compare(this->type()) == 0)
                  {
                     // Read version if present
                     getline(this->mFile, fileVers);

                     // Check that reading version was successful
                     if(this->mFile.good())
                     {
                        // Check version
                        if(fileVers.compare(this->version()) == 0)
                        {
                           //
                           // Compatibility check was successful
                           //
                        } else
                        {
                           throw Exception("Wrong ASCII file version!");
                        }
                     } else
                     {
                        throw Exception("Missing ASCII file version!");
                     }
                  } else
                  {
                     throw Exception("Wrong ASCII file type!");
                  }
               } else
               {
                  throw Exception("Missing ASCII file type!");
               }
            } else
            {
               throw Exception("Wrong ASCII file header!");
            }
         } else
         {
            throw Exception("Missing ASCII file header!");
         }
      }
   }

   void IAsciiReader::finalize()
   {
      // Check if the workflow allows IO to be performed
      if(FrameworkMacro::allowsIO())
      {
         // Close file
         this->close();
      }
   }

   void IAsciiReader::close()
   {
      // Check if the workflow allows IO to be performed
      if(FrameworkMacro::allowsIO())
      {
         this->mFile.close();
      }
   }

}
}
