/** 
 * @file BoussinesqRRB1DBoxVCContinuity.cpp
 * @brief Source of the implementation of the continuity equation in rotating Rayleigh-Benard convection in 1D box
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
#include "Equations/Box/Boussinesq/BoussinesqRRB1DBoxVCContinuity.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "Enums/NonDimensional.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   BoussinesqRRB1DBoxVCContinuity::BoussinesqRRB1DBoxVCContinuity(SharedEquationParameters spEqParams)
      : IScalarEquation(spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   BoussinesqRRB1DBoxVCContinuity::~BoussinesqRRB1DBoxVCContinuity()
   {
   }

   void BoussinesqRRB1DBoxVCContinuity::setCoupling()
   {
      this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::PROGNOSTIC, 1, false, true, false);
   }

   void BoussinesqRRB1DBoxVCContinuity::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const
   {
      // Assert on scalar component is used
      assert(id == FieldComponents::Physical::SCALAR);
   }

   void BoussinesqRRB1DBoxVCContinuity::setRequirements()
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
