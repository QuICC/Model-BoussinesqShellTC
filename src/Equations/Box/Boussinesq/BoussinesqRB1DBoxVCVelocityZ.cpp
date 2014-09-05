/** 
 * @file BoussinesqRB1DBoxVCVelocityZ.cpp
 * @brief Source of the implementation of th momentum equation for the Z component in Rayleigh-Benard convection in 1D box
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
#include "Equations/Box/Boussinesq/BoussinesqRB1DBoxVCVelocityZ.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "Enums/NonDimensional.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   BoussinesqRB1DBoxVCVelocityZ::BoussinesqRB1DBoxVCVelocityZ(SharedEquationParameters spEqParams)
      : IScalarEquation(spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   BoussinesqRB1DBoxVCVelocityZ::~BoussinesqRB1DBoxVCVelocityZ()
   {
   }

   void BoussinesqRB1DBoxVCVelocityZ::setCoupling()
   {
      this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::PROGNOSTIC, 1, false, true, false);
   }

   void BoussinesqRB1DBoxVCVelocityZ::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const
   {
      // Assert on scalar component is used
      assert(id == FieldComponents::Physical::SCALAR);
   }

   void BoussinesqRB1DBoxVCVelocityZ::setRequirements()
   {
      // Set temperatur as equation unknown
      this->setName(PhysicalNames::VELOCITYZ);

      // Set solver timing
      this->setSolveTiming(SolveTiming::PROGNOSTIC);

      // Add Z velocity to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::VELOCITYZ, FieldRequirement(true, true, true, true));

      // Add X velocity to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::VELOCITYX, FieldRequirement(true, true, true, false));

      // Add Y velocity to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::VELOCITYY, FieldRequirement(true, true, true, false));

      // Add temperature to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::TEMPERATURE, FieldRequirement(true, true, true, false));
   }

}
}
