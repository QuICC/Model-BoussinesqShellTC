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
      : mSteps(0), mAsciiRate(-1), mHdf5Rate(-1), mStatsRate(-1), mStatsAvgRate(-1), mActiveStatsUpdate(false), mActiveStatsWrite(true)
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

   bool SimulationIoControl::isAsciiTime() const
   {
      return (this->mAsciiRate > 0 && this->mSteps % this->mAsciiRate == 0);
   }

   bool SimulationIoControl::isHdf5Time() const
   {
      return (this->mHdf5Rate > 0 && this->mSteps % this->mHdf5Rate == 0);
   }

   void SimulationIoControl::activateStats()
   {
      bool writeTrigger = (this->mStatsRate > 0 && this->mSteps % this->mStatsRate == 0);
      bool updateTrigger = (this->mStatsAvgRate > 0 && this->mSteps % this->mStatsAvgRate == 0);

      this->mActiveStatsUpdate = (writeTrigger || updateTrigger);
   }

   void SimulationIoControl::disableStats()
   {
      this->mActiveStatsUpdate = false;
   }

   bool SimulationIoControl::isStatsTime() const
   {
      bool writeTrigger = (this->mStatsRate > 0 && this->mSteps % this->mStatsRate == 0);
      return writeTrigger;
   }

   bool SimulationIoControl::isStatsUpdateTime() const
   {
      return this->mActiveStatsUpdate;
   }

   void SimulationIoControl::writeFiles(const MHDFloat time, const MHDFloat timestep)
   {
      if(this->isAsciiTime())
      {
         this->writeAscii(time,timestep);
      }

      if(this->isHdf5Time())
      {
         this->writeHdf5(time,timestep);
      }

      // Activate stats
      this->activateStats();

      if(this->isStatsTime())
      {
         this->prepareStats(time,timestep);
         this->mActiveStatsWrite = true;
      }
   }

   void SimulationIoControl::finalize()
   {
      if(this->mspStdOut)
      {
         this->mspStdOut->finalize();
      }

      // Finalize ASCII and HDF5 writers
      this->finalizeWriters();
   }

   void SimulationIoControl::addAsciiOutputFile(IoVariable::SharedIVariableAsciiEWriter spOutFile)
   {
      this->mAsciiWriters.push_back(spOutFile);
   }

   void SimulationIoControl::addHdf5OutputFile(IoVariable::SharedIVariableHdf5NWriter spOutFile)
   {
      this->mHdf5Writers.push_back(spOutFile);
   }

   void SimulationIoControl::addStatsOutputFile(IoStats::SharedIStatisticsAsciiEWriter spOutFile)
   {
      this->mStatsWriters.push_back(spOutFile);
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

      // First check that all ASCII writers are full
      SimulationIoControl::stats_iterator itStats;
      for(itStats = this->mStatsWriters.begin(); itStats < this->mStatsWriters.end(); itStats++)
      {
         if(!(*itStats)->isFull())
         {
            throw Exception("There are missing variables in the statistics writers");
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

      // Iterate over all statistics writer
      for(itStats = this->mStatsWriters.begin(); itStats < this->mStatsWriters.end(); itStats++)
      {
         (*itStats)->setPhysical(phys);
         (*itStats)->init();
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
      this->mAsciiWriters.clear();

      // Iterate over all HDF5 writer
      SimulationIoControl::hdf5_iterator itHdf5;
      for(itHdf5 = this->mHdf5Writers.begin(); itHdf5 < this->mHdf5Writers.end(); itHdf5++)
      {
         (*itHdf5)->finalize();
      }
      this->mHdf5Writers.clear();

      // Iterate over all statistics writer
      SimulationIoControl::stats_iterator itStats;
      for(itStats = this->mStatsWriters.begin(); itStats < this->mStatsWriters.end(); itStats++)
      {
         (*itStats)->finalize();
      }
      this->mStatsWriters.clear();
   }

   void SimulationIoControl::writeAscii(const MHDFloat time, const MHDFloat timestep)
   {
      // Iterate over all ASCII writer
      SimulationIoControl::ascii_iterator it;
      for(it = this->mAsciiWriters.begin(); it < this->mAsciiWriters.end(); it++)
      {
         (*it)->setSimTime(time,timestep);
         (*it)->write();
      }
   }

   void SimulationIoControl::writeHdf5(const MHDFloat time, const MHDFloat timestep)
   {
      // Iterate over all HDF5 writer
      SimulationIoControl::hdf5_iterator it;
      for(it = this->mHdf5Writers.begin(); it < this->mHdf5Writers.end(); it++)
      {
         (*it)->setSimTime(time,timestep);
         (*it)->write();
      }
   }

   void SimulationIoControl::prepareStats(const MHDFloat time, const MHDFloat timestep)
   {
      // Iterate over all statistics writer
      SimulationIoControl::stats_iterator it;
      for(it = this->mStatsWriters.begin(); it < this->mStatsWriters.end(); it++)
      {
         (*it)->setSimTime(time,timestep);
      }
   }

   void SimulationIoControl::writeStats()
   {
      if(this->mActiveStatsWrite)
      {
         // Iterate over all statistics writer
         SimulationIoControl::stats_iterator it;
         for(it = this->mStatsWriters.begin(); it < this->mStatsWriters.end(); it++)
         {
            (*it)->write();
         }

         this->mActiveStatsWrite = false;
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

      // Set statistics rate from config file
      this->mStatsRate = this->mspCfgFile->spStats()->fValue("output_rate");
      this->mStatsAvgRate = this->mspCfgFile->spStats()->fValue("time_avg_rate");
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

   std::map<std::string, MHDFloat>& SimulationIoControl::rConfigPhysical()
   {
      // Safety assert for non NULL pointer
      assert(this->mspCfgFile);

      return this->mspCfgFile->rspPhysical()->rFMap();
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

   SimulationIoControl::stats_iterator  SimulationIoControl::beginStats()
   {
      return this->mStatsWriters.begin();
   }

   SimulationIoControl::stats_iterator  SimulationIoControl::endStats()
   {
      return this->mStatsWriters.end();
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

      Array tstep(3);

      // Get timestepping time configuration
      tstep(0) = this->mspCfgFile->spTimestepping()->fValue("time");

      // Get timestepping timestep configuration
      tstep(1) = this->mspCfgFile->spTimestepping()->fValue("timestep");

      // Get timestepping error configuration
      tstep(2) = this->mspCfgFile->spTimestepping()->fValue("error");

      return tstep;
   }
}
