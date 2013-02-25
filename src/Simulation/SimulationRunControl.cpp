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
      : mStatus(Runtime::Status::GOON), mCtrlFile(), mSteps(0)
   {
   }

   SimulationRunControl::~SimulationRunControl()
   {
   }

   Runtime::Status::Id SimulationRunControl::status() const
   {
      return this->mStatus;
   }

   bool SimulationRunControl::doIO() const
   {
      return (this->mSteps % 100) == 0;
   }

   void SimulationRunControl::update()
   {
      // Not ideal but OK for the moment
      this->mCtrlFile.read();

      this->mStatus = static_cast<Runtime::Status::Id>(static_cast<int>(this->mStatus) + static_cast<int>(this->mCtrlFile.status()));

      // Increment timestep counter
      this->mSteps++;
   }

   void SimulationRunControl::checkFile()
   {
      // Read input from control file
      this->mCtrlFile.read();
   }
}
