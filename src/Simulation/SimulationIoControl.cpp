/** 
 * @file SimulationIoControl.cpp
 * @brief Source of the implementation of the IO control system
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
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
      : mSteps(0), mAsciiRate(-1), mHdf5Rate(-1)
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

   void SimulationIoControl::update()
   {
      // Increment timestep counter
      this->mSteps++;
   }

   void SimulationIoControl::writeFiles(const MHDFloat time, const MHDFloat timestep)
   {
      if(this->mAsciiRate > 0 && this->mSteps % this->mAsciiRate == 0)
      {
         this->writeAscii(time,timestep);
      }

      if(this->mHdf5Rate > 0 && this->mSteps % this->mHdf5Rate == 0)
      {
         this->writeHdf5(time,timestep);
      }
   }

   void SimulationIoControl::finalize()
   {
      if(this->mspStdOut)
      {
         this->mspStdOut->finalize();
      }
   }

   void SimulationIoControl::addAsciiOutputFile(IoVariable::SharedIVariableAsciiEWriter spOutFile)
   {
      this->mAsciiWriters.push_back(spOutFile);
   }

   void SimulationIoControl::addHdf5OutputFile(IoVariable::SharedIVariableHdf5NWriter spOutFile)
   {
      this->mHdf5Writers.push_back(spOutFile);
   }

   void SimulationIoControl::initWriters()
   {  
      // First check that all ASCII writers are full
      SimulationIoControl::ascii_iterator itAscii;
      for(itAscii = this->mAsciiWriters.begin(); itAscii < this->mAsciiWriters.end(); itAscii++)
      {
         if(!(*itAscii)->isFull())
         {
            throw Exception("There are missing variables in the ASCII writers");
         }
      }

      // First check that all HDF5 writers are full
      SimulationIoControl::hdf5_iterator itHdf5;
      for(itHdf5 = this->mHdf5Writers.begin(); itHdf5 < this->mHdf5Writers.end(); itHdf5++)
      {
         if(!(*itHdf5)->isFull())
         {
            throw Exception("There are missing variables in the HDF5 writers");
         }
      }

      // Create physical parameters map
      std::map<std::string,MHDFloat> phys = this->configPhysical();
      std::map<std::string,int>::const_iterator it;
      for(it = this->configBoundary().begin(); it != this->configBoundary().end(); ++it)
      {
         phys.insert(std::make_pair("bc_"+it->first, it->second));
      }

      // Iterate over all ASCII writer
      for(itAscii = this->mAsciiWriters.begin(); itAscii < this->mAsciiWriters.end(); itAscii++)
      {
         (*itAscii)->setPhysical(phys);
         (*itAscii)->init();
      }

      // Iterate over all HDF5 writer
      for(itHdf5 = this->mHdf5Writers.begin(); itHdf5 < this->mHdf5Writers.end(); itHdf5++)
      {
         (*itHdf5)->setPhysical(phys);
         (*itHdf5)->init();
      }
   }

   void SimulationIoControl::finalizeWriters()
   {
      // Iterate over all ASCII writer
      SimulationIoControl::ascii_iterator itAscii;
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

   void SimulationIoControl::writeAscii(const MHDFloat time, const MHDFloat timestep)
   {
      // Iterate over all ASCII writer
      SimulationIoControl::ascii_iterator itAscii;
      for(itAscii = this->mAsciiWriters.begin(); itAscii < this->mAsciiWriters.end(); itAscii++)
      {
         (*itAscii)->setSimTime(time,timestep);
         (*itAscii)->write();
      }
   }

   void SimulationIoControl::writeHdf5(const MHDFloat time, const MHDFloat timestep)
   {
      // Iterate over all HDF5 writer
      SimulationIoControl::hdf5_iterator itHdf5;
      for(itHdf5 = this->mHdf5Writers.begin(); itHdf5 < this->mHdf5Writers.end(); itHdf5++)
      {
         (*itHdf5)->setSimTime(time,timestep);
         (*itHdf5)->write();
      }
   }

   void SimulationIoControl::initCfg()
   {
      // Safety assert for non NULL pointer
      assert(this->mspCfgFile);

      // Initialise config file
      this->mspCfgFile->init();

      // Read configuration data from config file
      this->mspCfgFile->read();

      // Close the config file
      this->mspCfgFile->finalize();

      // Set ASCII rate from config file
      this->mAsciiRate = this->mspCfgFile->spIo()->fValue("ascii");

      // Set HDF5 rate from config file
      this->mHdf5Rate = this->mspCfgFile->spIo()->fValue("hdf5");
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
      // Safety assert for non NULL pointer
      assert(this->mspCfgFile);

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

   Array SimulationIoControl::configBoxScale() const
   {
      // Safety assert for non NULL pointer
      assert(this->mspCfgFile);

      // Create storage for box scales
      Array box = Array::Ones(this->mspCfgFile->spTruncation()->iMap().size()); 

      // Get box scale truncation map
      std::map<std::string,MHDFloat>  trunc = this->mspCfgFile->spTruncation()->fMap();

      // Extract box scale from truncation read from file
      std::map<std::string,MHDFloat>::const_iterator  it;
      for(it = trunc.begin(); it != trunc.end(); it++)
      {
         if(it->first == "kc1D")
         {
            box(0) = it->second/trunc.find("box1D")->second;
         } else if(it->first == "kc2D")
         {
            box(1) = it->second/trunc.find("box2D")->second;
         } else if(it->first == "kc3D")
         {
            box(2) = it->second/trunc.find("box3D")->second;
         }
      }

      return box;
   }

   int SimulationIoControl::configNCpu() const
   {
      // Safety assert for non NULL pointer
      assert(this->mspCfgFile);

      return this->mspCfgFile->spParallel()->iValue("cpus");
   }

   const std::map<std::string, MHDFloat>& SimulationIoControl::configPhysical() const
   {
      // Safety assert for non NULL pointer
      assert(this->mspCfgFile);

      return this->mspCfgFile->spPhysical()->fMap();
   }

   const std::map<std::string, int>& SimulationIoControl::configBoundary() const
   {
      // Safety assert for non NULL pointer
      assert(this->mspCfgFile);

      return this->mspCfgFile->spBoundary()->iMap();
   }

   SimulationIoControl::ascii_iterator  SimulationIoControl::beginAscii()
   {
      return this->mAsciiWriters.begin();
   }

   SimulationIoControl::ascii_iterator  SimulationIoControl::endAscii()
   {
      return this->mAsciiWriters.end();
   }

   SimulationIoControl::hdf5_iterator  SimulationIoControl::beginHdf5()
   {
      return this->mHdf5Writers.begin();
   }

   SimulationIoControl::hdf5_iterator  SimulationIoControl::endHdf5()
   {
      return this->mHdf5Writers.end();
   }

   Array SimulationIoControl::configRun() const
   {
      // Safety assert for non NULL pointer
      assert(this->mspCfgFile);

      Array run(2);

      // Get simulation time configuration
      run(0) = this->mspCfgFile->spRun()->fValue("sim");

      // Get wall time configuration
      run(1) = this->mspCfgFile->spRun()->fValue("wall");

      return run;
   }

   Array SimulationIoControl::configTimestepping() const
   {
      // Safety assert for non NULL pointer
      assert(this->mspCfgFile);

      Array tstep(2);

      // Get timestepping time configuration
      tstep(0) = this->mspCfgFile->spTimestepping()->fValue("time");

      // Get timestepping timestep configuration
      tstep(1) = this->mspCfgFile->spTimestepping()->fValue("timestep");

      return tstep;
   }
}
