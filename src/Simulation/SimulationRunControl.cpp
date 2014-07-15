/** 
 * @file SimulationRunControl.cpp
 * @brief Implementation of a general simulation control structure
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// System includes
//
#include <csignal>

// External includes
//

// Class include
//
#include "Simulation/SimulationRunControl.hpp"

// Project includes
//
#include "IoTools/Formatter.hpp"

namespace GeoMHDiSCC {

   RuntimeStatus::Id SimulationRunControl::SIGNAL_STATUS = RuntimeStatus::GOON;

   SimulationRunControl::SimulationRunControl()
      : mStatus(RuntimeStatus::GOON), mCtrlFile(), mSteps(0), mMaxSimTime(0.0), mMaxWallTime(0.0)
   {
      // Initialise signal handling
      this->initSignalHandler();
   }

   SimulationRunControl::~SimulationRunControl()
   {
   }

   RuntimeStatus::Id SimulationRunControl::status() const
   {
      return this->mStatus;
   }

   void SimulationRunControl::updateSimulation(const MHDFloat simTime, const MHDFloat simDt)
   {
      // Increment timestep counter
      this->mSteps++;

      // Check for maximum simulation time
      if(this->mMaxSimTime > 0 && simTime > this->mMaxSimTime)
      {
         this->mStatus = RuntimeStatus::STOP;

         // Produce a nice looking output to std output 
         IoTools::Formatter::printLine(std::cout, '#');
         IoTools::Formatter::printCentered(std::cout, "Simulation time limit reached!", '#');
         IoTools::Formatter::printLine(std::cout, '#');
         IoTools::Formatter::printNewline(std::cout);
      }

      // Check for maximum simulation steps
      if(this->mMaxSimTime < 0 && this->mSteps > std::abs(this->mMaxSimTime))
      {
         this->mStatus = RuntimeStatus::STOP;

         // Produce a nice looking output to std output 
         IoTools::Formatter::printLine(std::cout, '#');
         IoTools::Formatter::printCentered(std::cout, "Simulation steps limit reached!", '#');
         IoTools::Formatter::printLine(std::cout, '#');
         IoTools::Formatter::printNewline(std::cout);
      }

      // Check if timestepper requested abort due to too small timestep
      if(simDt < 0)
      {
         this->mStatus = RuntimeStatus::STOP;

         // Produce a nice looking output to std output 
         IoTools::Formatter::printLine(std::cout, '#');
         IoTools::Formatter::printCentered(std::cout, "Adaptive timestep failed!", '#');
         IoTools::Formatter::printLine(std::cout, '#');
         IoTools::Formatter::printNewline(std::cout);
      }
   }

   void SimulationRunControl::updateCluster(const MHDFloat wallTime)
   {
      // Not ideal but OK for the moment
      this->mCtrlFile.read();

      // Convert wallTime to hours
      MHDFloat hours = wallTime/3600.;

      this->mStatus = static_cast<RuntimeStatus::Id>(static_cast<int>(this->mStatus) + static_cast<int>(this->mCtrlFile.status()));

      // Check for maximum wall time
      if(this->mMaxWallTime > 0 && hours > this->mMaxWallTime)
      {
         this->mStatus = RuntimeStatus::STOP;

         // Produce a nice looking output to std output 
         IoTools::Formatter::printLine(std::cout, '#');
         IoTools::Formatter::printCentered(std::cout, "Simulation wall time reached!", '#');
         IoTools::Formatter::printLine(std::cout, '#');
         IoTools::Formatter::printNewline(std::cout);
      }

      // Signal status
      this->mStatus = static_cast<RuntimeStatus::Id>(static_cast<int>(this->mStatus) + static_cast<int>(SimulationRunControl::SIGNAL_STATUS));
   }

   void SimulationRunControl::checkFile()
   {
      // Read input from control file
      this->mCtrlFile.read();
   }

   void SimulationRunControl::setMaxSimTime(const MHDFloat maxTime)
   {
      this->mMaxSimTime = maxTime;
   }

   void SimulationRunControl::setMaxWallTime(const MHDFloat maxTime)
   {
      this->mMaxWallTime = maxTime;
   }

   void SimulationRunControl::printInfo(std::ostream& stream)
   {
      // Create nice looking ouput header
      IoTools::Formatter::printNewline(stream);
      IoTools::Formatter::printLine(stream, '-');
      IoTools::Formatter::printCentered(stream, "Simulation run information", '*');
      IoTools::Formatter::printLine(stream, '-');

      std::stringstream oss;

      // get a nice base for info
      int base = 20;
      oss << std::fixed << std::setprecision(1) << this->mSteps;
      base += oss.str().size() + 1;
      oss.str("");

      // Output number of timesteps
      oss << "Timesteps: " << std::fixed << std::setprecision(1) << this->mSteps;
      IoTools::Formatter::printCentered(stream, oss.str(), ' ', base);
      oss.str("");

      IoTools::Formatter::printLine(stream, '*');
   }

   void SimulationRunControl::handleSignal(int signum)
   {
      if(signum == SIGUSR1)
      {
         SimulationRunControl::SIGNAL_STATUS = RuntimeStatus::STOP;

         IoTools::Formatter::printLine(std::cout, '#');
         IoTools::Formatter::printCentered(std::cout, "Simulation received stop signal!", '#');
         IoTools::Formatter::printLine(std::cout, '#');
      }
   }

   void SimulationRunControl::initSignalHandler()
   {
      signal(SIGUSR1, SimulationRunControl::handleSignal);
   }

}
