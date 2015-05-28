/** 
 * @file BoussinesqRBCCylinderVCContinuity.cpp
 * @brief Source of the implementation of the continuity equation in Rayleigh-Benard convection in a cylinder
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
#include "Equations/Cylinder/Boussinesq/BoussinesqRBCCylinderVCContinuity.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "Enums/NonDimensional.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   BoussinesqRBCCylinderVCContinuity::BoussinesqRBCCylinderVCContinuity(SharedEquationParameters spEqParams)
      : IScalarEquation(spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   BoussinesqRBCCylinderVCContinuity::~BoussinesqRBCCylinderVCContinuity()
   {
   }

   void BoussinesqRBCCylinderVCContinuity::setCoupling()
   {
      this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::PROGNOSTIC, 0, false, false);
   }

   void BoussinesqRBCCylinderVCContinuity::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const
   {
      // Assert on scalar component is used
      assert(id == FieldComponents::Physical::SCALAR);
   }

   void BoussinesqRBCCylinderVCContinuity::setRequirements()
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
