/** \file SimulationRunControl.cpp
 *  \brief Implementation of a general simulation control structure
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
      : mStatus(Runtime::Status::GOON), mCtrlFile(), mSteps(0), mMaxSimTime(0.0), mMaxWallTime(0.0)
   {
   }

   SimulationRunControl::~SimulationRunControl()
   {
   }

   Runtime::Status::Id SimulationRunControl::status() const
   {
      return this->mStatus;
   }

   void SimulationRunControl::update(const MHDFloat simTime)
   {
      // Increment timestep counter
      this->mSteps++;

      // Not ideal but OK for the moment
      this->mCtrlFile.read();

      this->mStatus = static_cast<Runtime::Status::Id>(static_cast<int>(this->mStatus) + static_cast<int>(this->mCtrlFile.status()));

      // Check for maximum simulation time
      if(this->mMaxSimTime > 0 && simTime > this->mMaxSimTime)
      {
         this->mStatus = Runtime::Status::STOP;
      }

      // Check for maximum simulation steps
      if(this->mMaxSimTime < 0 && this->mSteps > std::abs(this->mMaxSimTime))
      {
         this->mStatus = Runtime::Status::STOP;
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
