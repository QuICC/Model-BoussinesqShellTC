/** 
 * @file SimulationBoundary.cpp
 * @brief Implementation of a general simulation control structure
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "Simulation/SimulationBoundary.hpp"

// Project includes
//

namespace GeoMHDiSCC {

   SimulationBoundary::SimulationBoundary()
   {
   }

   SimulationBoundary::~SimulationBoundary()
   {
   }

   bool SimulationBoundary::hasEquation(const SpectralFieldId eqId) const
   {
      return this->mBcs.count(eqId);
   }

   bool SimulationBoundary::hasField(const SpectralFieldId eqId, const SpectralFieldId fieldId) const
   {
      // Safey assert
      assert(this->mBcs.count(eqId) > 0);

      return this->mBcs.find(eqId)->second.count(fieldId);
   }

   const SimulationBoundary::DimBCVectorMap& SimulationBoundary::bcs(const SpectralFieldId eqId, const SpectralFieldId fieldId) const
   {
      // Safey asserts
      assert(this->mBcs.count(eqId) > 0);
      assert(this->mBcs.find(eqId)->second.count(fieldId) > 0);

      return this->mBcs.find(eqId)->second.find(fieldId)->second;
   }

   void SimulationBoundary::initStorage(const SpectralFieldId id)
   {
      // Initialise boundary condition storage
      this->mBcs.insert(std::make_pair(id, std::map<SpectralFieldId, SimulationBoundary::DimBCVectorMap>()));
   }

   void SimulationBoundary::initBcStorage(const SpectralFieldId eqId, const SpectralFieldId fieldId, const Dimensions::Simulation::Id dimId)
   {
      // Safey assert
      assert(this->mBcs.count(eqId) > 0);

      // Initialise boundary condition field storage
      this->mBcs.find(eqId)->second.insert(std::make_pair(fieldId, SimulationBoundary::DimBCVectorMap()));

      // Initialise boundary condition storage
      this->mBcs.find(eqId)->second.find(fieldId)->second.insert(std::make_pair(dimId, Boundary::BCVector()));
   }

   void SimulationBoundary::addBc(const SpectralFieldId eqId, const SpectralFieldId fieldId, const Dimensions::Simulation::Id dimId, Boundary::BCType bcId, Boundary::BCPosition pos)
   {
      // Safey assert
      assert(this->mBcs.count(eqId) > 0);
      assert(this->mBcs.find(eqId)->second.count(fieldId) > 0);
      assert(this->mBcs.find(eqId)->second.find(fieldId)->second.count(dimId) > 0);

      this->mBcs.find(eqId)->second.find(fieldId)->second.find(dimId)->second.push_back(Boundary::BoundaryCondition(bcId,pos));
   }

}
