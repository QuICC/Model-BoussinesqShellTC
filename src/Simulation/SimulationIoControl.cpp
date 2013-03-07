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

   void SimulationIoControl::addOutputFile(IoVariable::SharedIVariableHdf5NWriter spOutFile)
   {
      this->mHdf5Writers.push_back(spOutFile);
   }

   void SimulationIoControl::initWriters()
   {  
      // First check that all HDF5 writers are full
      SimulationIoControl::hdf5_iterator itHdf5;
      for(itHdf5 = this->mHdf5Writers.begin(); itHdf5 < this->mHdf5Writers.end(); itHdf5++)
      {
         if(!(*itHdf5)->isFull())
         {
            throw Exception("There are missing variables in the HDF5 writers");
         }
      }

      // Iterate over all ASCII writer
      std::vector<IoAscii::SharedIAsciiWriter>::iterator itAscii;
      for(itAscii = this->mAsciiWriters.begin(); itAscii < this->mAsciiWriters.end(); itAscii++)
      {
         (*itAscii)->init();
      }

      // Iterate over all HDF5 writer
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
      SimulationIoControl::hdf5_iterator itHdf5;
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
      SimulationIoControl::hdf5_iterator itHdf5;
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

   ArrayI SimulationIoControl::configDimension() const
   {
      // Get truncation map
      std::map<std::string,int>  trunc = this->mspCfgFile->spTruncation()->iMap();

      // Create storage for the dimensions
      ArrayI dim(trunc.size());

      // Extrac dimension from truncation read from file
      std::map<std::string,int>::const_iterator  itI;
      int i = 0;
      for(itI = trunc.begin(); itI != trunc.end(); itI++)
      {
         dim(i) = itI->second;
         i++;
      }

      return dim;
   }

   int SimulationIoControl::configNCpu() const
   {
      return this->mspCfgFile->spParallel()->iValue("cpus");
   }

   const std::map<std::string, MHDFloat>& SimulationIoControl::configPhysical() const
   {
      return this->mspCfgFile->spPhysical()->fMap();
   }

   const std::map<std::string, int>& SimulationIoControl::configBoundary() const
   {
      return this->mspCfgFile->spBoundary()->iMap();
   }

   SimulationIoControl::hdf5_iterator  SimulationIoControl::beginHdf5()
   {
      return this->mHdf5Writers.begin();
   }

   SimulationIoControl::hdf5_iterator  SimulationIoControl::endHdf5()
   {
      return this->mHdf5Writers.end();
   }
}
