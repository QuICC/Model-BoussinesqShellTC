/** \file RandomScalarState.cpp
 *  \brief Source of the implementation of the general random scalar state equation
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "Generator/States/RandomScalarState.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "TypeSelectors/SpectralSelector.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   RandomScalarState::RandomScalarState(SharedEquationParameters spEqParams, const PhysicalNames::Id name)
      : IScalarEquation(spEqParams)
   {
      // Set name of unknown
      this->setName(name);

      // Set the variable requirements
      this->setRequirements();
   }

   RandomScalarState::~RandomScalarState()
   {
   }

   void RandomScalarState::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const
   {
      // Assert on scalar component is used
      assert(id == FieldComponents::Physical::SCALAR);

   }

   void RandomScalarState::setRequirements()
   {
      // Add unknown to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(this->name(), FieldRequirement(true, true, false, false));
   }

   void RandomScalarState::setCoupling()
   {
   }

}
}
