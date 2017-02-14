/** 
 * @file BoussinesqRRBCBoxVCMomentum.cpp
 * @brief Source of the implementation of the vector momentum equation for rotating Rayleigh-Benard convection in 3D box
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
#include "Equations/Box/Boussinesq/BoussinesqRRBCBoxVCMomentum.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "Enums/NonDimensional.hpp"
#include "PhysicalOperators/VelocityAdvection.hpp"

namespace QuICC {

namespace Equations {

   BoussinesqRRBCBoxVCMomentum::BoussinesqRRBCBoxVCMomentum(SharedEquationParameters spEqParams)
      : IVectorEquation(spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   BoussinesqRRBCBoxVCMomentum::~BoussinesqRRBCBoxVCMomentum()
   {
   }

   void BoussinesqRRBCBoxVCMomentum::setCoupling()
   {
      this->defineCoupling(FieldComponents::Spectral::X, CouplingInformation::PROGNOSTIC, 0, true, false);

      this->defineCoupling(FieldComponents::Spectral::Y, CouplingInformation::PROGNOSTIC, 0, true, false);

      this->defineCoupling(FieldComponents::Spectral::Z, CouplingInformation::PROGNOSTIC, 0, true, false);
   }

   void BoussinesqRRBCBoxVCMomentum::setNLComponents()
   {
      this->addNLComponent(FieldComponents::Spectral::X, 0);

      this->addNLComponent(FieldComponents::Spectral::Y, 0);

      this->addNLComponent(FieldComponents::Spectral::Z, 0);
   }

   void BoussinesqRRBCBoxVCMomentum::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const
   {
      // Get Prandtl number
      MHDFloat Pr = this->eqParams().nd(NonDimensional::PRANDTL);

      /// 
      /// Computation of the advection:
      ///   \f$ \left(\vec u\cdot\nabla\right)u_i\f$
      ///
      FieldComponents::Spectral::Id specId = FieldComponents::Spectral::NOTUSED;
      if(id == FieldComponents::Physical::X)
      {
         specId = FieldComponents::Spectral::X;
      } else if(id == FieldComponents::Physical::Y)
      {
         specId = FieldComponents::Spectral::Y;
      } else if(id == FieldComponents::Physical::Z)
      {
         specId = FieldComponents::Spectral::Z;
      }

      Physical::VelocityAdvection<FieldComponents::Physical::X,FieldComponents::Physical::Y,FieldComponents::Physical::Z>::set(rNLComp, this->unknown().dom(0).phys(), this->unknown().dom(0).grad(specId), 1.0/Pr);
   }

   void BoussinesqRRBCBoxVCMomentum::setRequirements()
   {
      // Set temperatur as equation unknown
      this->setName(PhysicalNames::VELOCITY);

      // Set solver timing
      this->setSolveTiming(SolveTiming::PROGNOSTIC);

      // Add velocity to requirements: is scalar?, need spectral?, need physical?, need grad?(, need curl?)
      this->mRequirements.addField(PhysicalNames::VELOCITY, FieldRequirement(false, true, true, true));
   }

}
}
