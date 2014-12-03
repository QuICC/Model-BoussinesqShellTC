/** 
 * @file BoussinesqRB3DBoxVCVelocity.cpp
 * @brief Source of the implementation of the vector momentum equation in Rayleigh-Benard convection in 3D box
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
#include "Equations/Box/Boussinesq/BoussinesqRB3DBoxVCVelocity.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "Enums/NonDimensional.hpp"
#include "PhysicalOperators/VelocityAdvection.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   BoussinesqRB3DBoxVCVelocity::BoussinesqRB3DBoxVCVelocity(SharedEquationParameters spEqParams)
      : IVectorEquation(spEqParams)
   {
      // Vector equation has three components
      this->mSpectralIds.push_back(FieldComponents::Spectral::ONE);
      this->mSpectralIds.push_back(FieldComponents::Spectral::TWO);
      this->mSpectralIds.push_back(FieldComponents::Spectral::THREE);

      // Set the variable requirements
      this->setRequirements();
   }

   BoussinesqRB3DBoxVCVelocity::~BoussinesqRB3DBoxVCVelocity()
   {
   }

   void BoussinesqRB3DBoxVCVelocity::setCoupling()
   {
      this->defineCoupling(FieldComponents::Spectral::ONE, CouplingInformation::PROGNOSTIC, 0, true, true, false);

      this->defineCoupling(FieldComponents::Spectral::TWO, CouplingInformation::PROGNOSTIC, 0, true, true, false);

      this->defineCoupling(FieldComponents::Spectral::THREE, CouplingInformation::PROGNOSTIC, 0, true, true, false);
   }

   void BoussinesqRB3DBoxVCVelocity::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const
   {
      /// 
      /// Computation of the advection:
      ///   \f$ \left(\vec u\cdot\nabla\right)u_i\f$
      ///
      Physical::VelocityAdvection<FieldComponents::Physical::ONE,FieldComponents::Physical::TWO,FieldComponents::Physical::THREE>::set(rNLComp, this->unknown().dom(0).phys().comp(FieldComponents::Physical::ONE), this->unknown().dom(0).phys().comp(FieldComponents::Physical::TWO), this->unkwnown().dom(0).phys().comp(FieldComponents::Physical::ONE), this->unknown().dom(0).grad().comp(id), 1.0);
   }

   void BoussinesqRB3DBoxVCVelocity::setRequirements()
   {
      // Set temperatur as equation unknown
      this->setName(PhysicalNames::VELOCITY);

      // Set solver timing
      this->setSolveTiming(SolveTiming::PROGNOSTIC);

      // Add velocity to requirements: is scalar?, need spectral?, need physical?, need grad?
      this->mRequirements.addField(PhysicalNames::VELOCITYX, FieldRequirement(false, true, true, true));
   }

}
}
