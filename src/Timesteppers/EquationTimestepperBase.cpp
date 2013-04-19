/** \file EquationTimestepperBase.cpp
 *  \brief Implementation of the base for an equation timestepper
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "Timesteppers/EquationTimestepperBase.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Timestep {

   EquationTimestepperBase::EquationTimestepperBase(const int start)
      : mZeroIdx(start)
   {
      // Safety assert
      assert(start >= 0);
   }

   EquationTimestepperBase::~EquationTimestepperBase()
   {
   }

   void EquationTimestepperBase::addInformation(const std::pair<PhysicalNames::Id, FieldComponents::Spectral::Id>& id, const ArrayI& startRow)
   {
      this->mInformation.insert(std::make_pair(id, startRow));
   }

   int EquationTimestepperBase::startRow(const std::pair<PhysicalNames::Id, FieldComponents::Spectral::Id>& id, const int i) const
   {
      return this->mInformation.find(id)->second(i);
   }
}
}
