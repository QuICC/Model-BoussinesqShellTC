/** \file ControlInterface.cpp
 *  \brief Source of the external control interface implementation
 */

// Configuration includes
//
#include "Framework/FrameworkMacro.h"

// System includes
//
#include <iostream>
#include <sstream>

// External includes
//

// Class include
//
#include "IoControl/ControlInterface.hpp"

// Project includes
//
#include "Exceptions/Exception.hpp"
#include "IoTools/Formatter.hpp"

namespace GeoMHDiSCC {

namespace IoControl {


   ControlInterface::ControlInterface()
      : ControlFile(), mStatus(Runtime::Status::GOON)
   {
      this->init();
   }

   ControlInterface::~ControlInterface()
   {
      this->finalize();
   }

   Runtime::Status::Id  ControlInterface::status() const
   {
      return this->mStatus;
   }

   void ControlInterface::init()
   {
      // Check if the framework allows IO to be performed
      if(FrameworkMacro::allowsIO())
      {
         // Create file
         this->create();

         // Close file
         this->close();
      }
   }

   void ControlInterface::read()
   {
      // Reset the current status
      this->mStatus = Runtime::Status::GOON;

      // Check if the framework allows IO to be performed
      if(FrameworkMacro::allowsIO())
      {
         // Open file
         this->open();

         // Update current status
         this->update();

         // Close file
         this->close();
      }

      // Share the status over all MPI processes
      #ifdef GEOMHDISCC_MPI
         MPI_Bcast(&this->mStatus, 1, MPI_INT, FrameworkMacro::IO_RANK, MPI_COMM_WORLD);
      #endif // GEOMHDISCC_MPI
   }

   void ControlInterface::finalize()
   {
      // Check if the framework allows IO to be performed
      if(FrameworkMacro::allowsIO())
      {
         // Close the file
         this->close();

         // Delete the file
         this->deleteFile();
      }
   }

   void ControlInterface::create()
   {
      // Check if the framework allows IO to be performed
      if(FrameworkMacro::allowsIO())
      {
         // Create file
         this->mFile.open(this->filename().c_str(), std::fstream::out);

         // Check if creation was successful
         if(! this->mFile.is_open())
         {
            throw Exception("Couldn't create control file " + this->filename() + "!");
         }
      }
   }

   void ControlInterface::open()
   {
      // Open file
      this->mFile.open(filename().c_str());

      // Check if file could be opened
      if(! this->mFile.is_open())
      {
         // File doesn't exits abort run
         this->mStatus = Runtime::Status::STOP;

      // File exist keep going
      } else
      {
         this->mStatus = Runtime::Status::GOON;
      }
   }

   void ControlInterface::close()
   {
      // Check if the framework allows IO to be performed
      if(FrameworkMacro::allowsIO())
      {
         this->mFile.close();
      }
   }

   void ControlInterface::deleteFile()
   {
      // Delete the file
      std::remove(this->filename().c_str());
   }

   void ControlInterface::update()
   {
      // Check if status requests stop of simulation
      if(this->mStatus == Runtime::Status::STOP)
      {
         // Check if the framework allows IO to be performed
         if(FrameworkMacro::allowsIO())
         {
            // Produce a nice looking output to std output 
            IoTools::Formatter::printLine(std::cout, '#');
            IoTools::Formatter::printCentered(std::cout, "User requested abort!", '#');
            IoTools::Formatter::printLine(std::cout, '#');
         }
      }
   }

}
}