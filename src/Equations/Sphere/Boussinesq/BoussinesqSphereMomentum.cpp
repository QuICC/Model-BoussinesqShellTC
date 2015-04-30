/** 
 * @file BoussinesqSphereMomentum.cpp
 * @brief Source of the implementation of the vector Navier-Stokes equation in the Boussinesq sphere model
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
#include "Equations/Sphere/Boussinesq/BoussinesqSphereMomentum.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "TypeSelectors/EquationEigenSelector.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   BoussinesqSphereMomentum::BoussinesqSphereMomentum(SharedEquationParameters spEqParams)
      : IVectorEquation(spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   BoussinesqSphereMomentum::~BoussinesqSphereMomentum()
   {
   }

   void BoussinesqSphereMomentum::setCoupling()
   {
      this->defineCoupling(FieldComponents::Spectral::ONE, CouplingInformation::PROGNOSTIC, 1, false, false);

      this->defineCoupling(FieldComponents::Spectral::TWO, CouplingInformation::PROGNOSTIC, 1, false, false);
   }

   void BoussinesqSphereMomentum::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const
   {
      throw Exception("Nonlinear term in spherical shell toroidal/poloidal Navier-Stokes equation not done yet");
   }

   void BoussinesqSphereMomentum::setRequirements()
   {
      // Set velocity as equation unknown
      this->setName(PhysicalNames::VELOCITY);

      // Set solver timing
      this->setSolveTiming(SolveTiming::PROGNOSTIC);

      // Add velocity to requirements: is scalar?, need spectral?, need physical?, need diff?(,need curl?)
      this->mRequirements.addField(PhysicalNames::VELOCITY, FieldRequirement(false, true, false, false));
   }

}
}
