/** 
 * @file BoussinesqCylinderVelocity.cpp
 * @brief Source of the implementation of the Navier-Stokes equation in the Boussinesq cylinder model
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
#include "Equations/Cylinder/Boussinesq/BoussinesqCylinderVelocity.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "TypeSelectors/SpectralOperatorSelector.hpp"
#include "TypeSelectors/EquationEigenSelector.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   BoussinesqCylinderVelocity::BoussinesqCylinderVelocity(SharedEquationParameters spEqParams)
      : IVectorEquation(spEqParams)
   {
      // Vector equation has two components, ie Toroidal/Poloidal
      this->mSpectralIds.push_back(FieldComponents::Spectral::ONE);
      this->mSpectralIds.push_back(FieldComponents::Spectral::TWO);

      // Set the variable requirements
      this->setRequirements();
   }

   BoussinesqCylinderVelocity::~BoussinesqCylinderVelocity()
   {
   }

   void BoussinesqCylinderVelocity::setCoupling()
   {
      this->defineCoupling(FieldComponents::Spectral::ONE, CouplingInformation::PROGNOSTIC, 0, false, false, false);

      this->defineCoupling(FieldComponents::Spectral::TWO, CouplingInformation::PROGNOSTIC, 0, false, false, false);
   }

   void BoussinesqCylinderVelocity::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const
   {
      throw Exception("Nonlinear term in spherical shell toroidal/poloidal Navier-Stokes equation not done yet");
   }

   void BoussinesqCylinderVelocity::setRequirements()
   {
      // Set velocity as equation unknown
      this->setName(PhysicalNames::VELOCITY);

      // Add velocity to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::VELOCITY, FieldRequirement(false, true, false, false));
   }

}
}
