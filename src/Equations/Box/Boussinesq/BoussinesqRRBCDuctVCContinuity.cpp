/** 
 * @file BoussinesqRRBCDuctVCContinuity.cpp
 * @brief Source of the implementation of the continuity equation in rotating Rayleigh-Benard convection in a infinite duct
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
#include "Equations/Box/Boussinesq/BoussinesqRRBCDuctVCContinuity.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "Enums/NonDimensional.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   BoussinesqRRBCDuctVCContinuity::BoussinesqRRBCDuctVCContinuity(SharedEquationParameters spEqParams)
      : IScalarEquation(spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   BoussinesqRRBCDuctVCContinuity::~BoussinesqRRBCDuctVCContinuity()
   {
   }

   void BoussinesqRRBCDuctVCContinuity::setCoupling()
   {
      this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::PROGNOSTIC, 0, false, false);
   }

   void BoussinesqRRBCDuctVCContinuity::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const
   {
      // Assert on scalar component is used
      assert(id == FieldComponents::Physical::SCALAR);
   }

   void BoussinesqRRBCDuctVCContinuity::setRequirements()
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
