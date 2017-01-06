/** 
 * @file BoussinesqRBCAnnulusVCContinuity.cpp
 * @brief Source of the implementation of the continuity equation in Rayleigh-Benard convection in a cylindrical annulus
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
#include "Equations/Annulus/Boussinesq/BoussinesqRBCAnnulusVCContinuity.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "Enums/NonDimensional.hpp"

namespace QuICC {

namespace Equations {

   BoussinesqRBCAnnulusVCContinuity::BoussinesqRBCAnnulusVCContinuity(SharedEquationParameters spEqParams)
      : IScalarEquation(spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   BoussinesqRBCAnnulusVCContinuity::~BoussinesqRBCAnnulusVCContinuity()
   {
   }

   void BoussinesqRBCAnnulusVCContinuity::setCoupling()
   {
      this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::PROGNOSTIC, 0, false, false);
   }

   void BoussinesqRBCAnnulusVCContinuity::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const
   {
      // Assert on scalar component is used
      assert(id == FieldComponents::Physical::SCALAR);
   }

   void BoussinesqRBCAnnulusVCContinuity::setRequirements()
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
