/** 
 * @file BoussinesqRB2DBoxVCContinuity.cpp
 * @brief Source of the implementation of the continuity equation in Rayleigh-Benard convection in 2D box
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "Equations/Box/Boussinesq/BoussinesqRB2DBoxVCContinuity.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "Enums/NonDimensional.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   BoussinesqRB2DBoxVCContinuity::BoussinesqRB2DBoxVCContinuity(SharedEquationParameters spEqParams)
      : IScalarEquation(spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   BoussinesqRB2DBoxVCContinuity::~BoussinesqRB2DBoxVCContinuity()
   {
   }

   void BoussinesqRB2DBoxVCContinuity::setCoupling()
   {
      this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::PROGNOSTIC, 0, false, true, false);
   }

   void BoussinesqRB2DBoxVCContinuity::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const
   {
      // Assert on scalar component is used
      assert(id == FieldComponents::Physical::SCALAR);
   }

   void BoussinesqRB2DBoxVCContinuity::setRequirements()
   {
      // Set temperatur as equation unknown
      this->setName(PhysicalNames::PRESSURE);

      // Set solver timing
      this->setSolveTiming(SolveTiming::PROGNOSTIC);

      // Add pressure to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::PRESSURE, FieldRequirement(true, true, false, false));
   }

}
}
