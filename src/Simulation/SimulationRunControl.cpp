/** 
 * @file SimulationRunControl.cpp
 * @brief Implementation of a general simulation control structure
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// System includes
//

// External includes
//

// Class include
//
#include "Simulation/SimulationRunControl.hpp"

// Project includes
//
#include "IoTools/Formatter.hpp"

namespace GeoMHDiSCC {

   SimulationRunControl::SimulationRunControl()
      : mStatus(RuntimeStatus::GOON), mCtrlFile(), mSteps(0), mMaxSimTime(0.0), mMaxWallTime(0.0)
   {
   }

   SimulationRunControl::~SimulationRunControl()
   {
   }

   RuntimeStatus::Id SimulationRunControl::status() const
   {
      return this->mStatus;
   }

   void SimulationRunControl::update(const MHDFloat simTime, const MHDFloat simDt)
   {
      // Increment timestep counter
      this->mSteps++;

      // Not ideal but OK for the moment
      this->mCtrlFile.read();

      this->mStatus = static_cast<RuntimeStatus::Id>(static_cast<int>(this->mStatus) + static_cast<int>(this->mCtrlFile.status()));

      // Check for maximum simulation time
      if(this->mMaxSimTime > 0 && simTime > this->mMaxSimTime)
      {
         this->mStatus = RuntimeStatus::STOP;

         // Produce a nice looking output to std output 
         IoTools::Formatter::printLine(std::cout, '#');
         IoTools::Formatter::printCentered(std::cout, "Simulation time limit reached!", '#');
         IoTools::Formatter::printLine(std::cout, '#');
      }

      // Check for maximum simulation steps
      if(this->mMaxSimTime < 0 && this->mSteps > std::abs(this->mMaxSimTime))
      {
         this->mStatus = RuntimeStatus::STOP;

         // Produce a nice looking output to std output 
         IoTools::Formatter::printLine(std::cout, '#');
         IoTools::Formatter::printCentered(std::cout, "Simulation steps limit reached!", '#');
         IoTools::Formatter::printLine(std::cout, '#');
      }

      // Check if timestepper requested abort due to too small timestep
      if(simDt < 0)
      {
         this->mStatus = RuntimeStatus::STOP;

         // Produce a nice looking output to std output 
         IoTools::Formatter::printLine(std::cout, '#');
         IoTools::Formatter::printCentered(std::cout, "Adaptive timestep failed!", '#');
         IoTools::Formatter::printLine(std::cout, '#');
      }
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
}
