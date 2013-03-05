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

   const Equations::IEvolutionEquation::BcEqMapType& SimulationBoundary::bc(const PhysicalNames::Id id) const
   {
      // Safey assert
      assert(this->mBcs.count(id) > 0);

      return this->mBcs.find(id)->second;
   }

   const std::map<PhysicalNames::Id, Equations::IEvolutionEquation::BcEqMapType>& SimulationBoundary::cbc(const PhysicalNames::Id id) const
   {
      // Safey assert
      assert(this->mCBcs.count(id) > 0);

      return this->mCBcs.find(id)->second;
   }

   void SimulationBoundary::initStorage(const PhysicalNames::Id id)
   {
      // Initialise boundary condition storage
      this->mBcs.insert(std::make_pair(id, Equations::IEvolutionEquation::BcEqMapType()));

      // Initialise coupled boundary condition storage
      this->mCBcs.insert(std::make_pair(id, std::map<PhysicalNames::Id, Equations::IEvolutionEquation::BcEqMapType>()));
   }

   void SimulationBoundary::initBcStorage(const PhysicalNames::Id id, const Equations::IEvolutionEquation::BcKeyType& key)
   {
      // Safey assert
      assert(this->mBcs.count(id) > 0);

      // Initialise boundary condition storage
      this->mBcs.find(id)->second.insert(std::make_pair(key, Equations::IEvolutionEquation::BcMapType()));
   }

   void SimulationBoundary::initCBcStorage(const PhysicalNames::Id eqId, const PhysicalNames::Id cpId, const Equations::IEvolutionEquation::BcKeyType& key)
   {
      // Safey assert
      assert(this->mCBcs.count(eqId) > 0);

      // Initialise coupled boundary condition field storage
      this->mCBcs.find(eqId)->second.insert(std::make_pair(cpId, Equations::IEvolutionEquation::BcEqMapType()));

      // Initialise boundary condition storage
      this->mCBcs.find(eqId)->second.find(cpId)->second.insert(std::make_pair(key, Equations::IEvolutionEquation::BcMapType()));
   }

   void SimulationBoundary::addBc(const PhysicalNames::Id, const Equations::IEvolutionEquation::BcKeyType& key, Spectral::BoundaryConditions::Id bcId, Spectral::IBoundary::Position pos)
   {
   }

   void SimulationBoundary::addCBc(const PhysicalNames::Id eqId, const PhysicalNames::Id cpId, const Equations::IEvolutionEquation::BcKeyType& key, Spectral::BoundaryConditions::Id bcId, Spectral::IBoundary::Position pos)
   {
   }

}
