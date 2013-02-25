/** \file SimulationIoControl.cpp
 *  \brief Source of the implementation of the IO control system
 */

// System includes
//

// External includes
//

// Class include
//
#include "Simulation/SimulationIoControl.hpp"

// Project includes
//

namespace GeoMHDiSCC {

   SimulationIoControl::SimulationIoControl()
   {
   }

   SimulationIoControl::~SimulationIoControl()
   {
   }

   void SimulationIoControl::setConfigurationFile(IoConfig::SharedConfigurationReader spCfgFile)
   {
      this->mspCfgFile = spCfgFile;
   }

   void SimulationIoControl::init()
   {
      // Init configuration file
      this->initCfg();

      // Init StdOutPipe output file
      this->initStdOut();

      // Print parameters file parameters
      this->mspCfgFile->printInfo();
   }

   void SimulationIoControl::cleanup()
   {
      // reset the configuration file
      this->mspCfgFile.reset();
   }

   void SimulationIoControl::finalize()
   {
      if(this->mspStdOut)
      {
         this->mspStdOut->finalize();
      }
   }

   void SimulationIoControl::addOutputFile(IoAscii::SharedIAsciiWriter spOutFile)
   {
      this->mAsciiWriters.push_back(spOutFile);
   }

   void SimulationIoControl::addOutputFile(IoHdf5::SharedIHdf5Writer spOutFile)
   {
      this->mHdf5Writers.push_back(spOutFile);
   }

   void SimulationIoControl::initWriters()
   {
      // Iterate over all ASCII writer
      std::vector<IoAscii::SharedIAsciiWriter>::iterator itAscii;
      for(itAscii = this->mAsciiWriters.begin(); itAscii < this->mAsciiWriters.end(); itAscii++)
      {
         (*itAscii)->init();
      }

      // Iterate over all HDF5 writer
      std::vector<IoHdf5::SharedIHdf5Writer>::iterator itHdf5;
      for(itHdf5 = this->mHdf5Writers.begin(); itHdf5 < this->mHdf5Writers.end(); itHdf5++)
      {
         (*itHdf5)->init();
      }
   }

   void SimulationIoControl::finalizeWriters()
   {
      // Iterate over all ASCII writer
      std::vector<IoAscii::SharedIAsciiWriter>::iterator itAscii;
      for(itAscii = this->mAsciiWriters.begin(); itAscii < this->mAsciiWriters.end(); itAscii++)
      {
         (*itAscii)->finalize();
      }

      // Iterate over all HDF5 writer
      std::vector<IoHdf5::SharedIHdf5Writer>::iterator itHdf5;
      for(itHdf5 = this->mHdf5Writers.begin(); itHdf5 < this->mHdf5Writers.end(); itHdf5++)
      {
         (*itHdf5)->finalize();
      }
   }

   void SimulationIoControl::writeAscii()
   {
      // Iterate over all ASCII writer
      std::vector<IoAscii::SharedIAsciiWriter>::iterator itAscii;
      for(itAscii = this->mAsciiWriters.begin(); itAscii < this->mAsciiWriters.end(); itAscii++)
      {
         (*itAscii)->write();
      }
   }

   void SimulationIoControl::writeHdf5()
   {
      // Iterate over all HDF5 writer
      std::vector<IoHdf5::SharedIHdf5Writer>::iterator itHdf5;
      for(itHdf5 = this->mHdf5Writers.begin(); itHdf5 < this->mHdf5Writers.end(); itHdf5++)
      {
         (*itHdf5)->write();
      }
   }

   void SimulationIoControl::initCfg()
   {
      // Initialise config file
      this->mspCfgFile->init();

      // Read configuration data from config file
      this->mspCfgFile->read();

      // Close the config file
      this->mspCfgFile->finalize();
   }

   void SimulationIoControl::initStdOut()
   {
      // Create the StdOutPipe ouput
      IoAscii::SharedStdOutPipe spStd(new IoAscii::StdOutPipe("OUT"));

      // Store the shared pointer
      this->mspStdOut = spStd;

      // Initialise the StdOutPipe output
      this->mspStdOut->init(); 
   }
}
