/** \file SimulationBoundary.cpp
 *  \brief Implementation of a general simulation control structure
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

   bool SimulationBoundary::hasEquation(const PhysicalNames::Id eqId) const
   {
      return this->mBcs.count(eqId);
   }

   bool SimulationBoundary::hasField(const PhysicalNames::Id eqId, const PhysicalNames::Id fieldId) const
   {
      // Safey assert
      assert(this->mBcs.count(eqId) > 0);

      return this->mBcs.find(eqId)->second.count(fieldId);
   }

   const SimulationBoundary::BcEqMapType& SimulationBoundary::bcs(const PhysicalNames::Id eqId, const PhysicalNames::Id fieldId) const
   {
      // Safey asserts
      assert(this->mBcs.count(eqId) > 0);
      assert(this->mBcs.find(eqId)->second.count(fieldId) > 0);

      return this->mBcs.find(eqId)->second.find(fieldId)->second;
   }

   void SimulationBoundary::initStorage(const PhysicalNames::Id id)
   {
      // Initialise boundary condition storage
      this->mBcs.insert(std::make_pair(id, std::map<PhysicalNames::Id, SimulationBoundary::BcEqMapType>()));
   }

   void SimulationBoundary::initBcStorage(const PhysicalNames::Id eqId, const PhysicalNames::Id fieldId, const SimulationBoundary::BcKeyType& key)
   {
      // Safey assert
      assert(this->mBcs.count(eqId) > 0);

      // Initialise boundary condition field storage
      this->mBcs.find(eqId)->second.insert(std::make_pair(fieldId, SimulationBoundary::BcEqMapType()));

      // Initialise boundary condition storage
      this->mBcs.find(eqId)->second.find(fieldId)->second.insert(std::make_pair(key, SimulationBoundary::BcMapType()));
   }

   void SimulationBoundary::addBc(const PhysicalNames::Id eqId, const PhysicalNames::Id fieldId, const SimulationBoundary::BcKeyType& key, Spectral::BoundaryConditions::Id bcId, Spectral::IBoundary::Position pos)
   {
      // Safey assert
      assert(this->mBcs.count(eqId) > 0);
      assert(this->mBcs.find(eqId)->second.count(fieldId) > 0);
      assert(this->mBcs.find(eqId)->second.find(fieldId)->second.count(key) > 0);

      this->mBcs.find(eqId)->second.find(fieldId)->second.find(key)->second.push_back(std::make_pair(bcId,pos));
   }

}
