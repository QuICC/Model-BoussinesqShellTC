/** 
 * @file BoussinesqRRBCAnnulusVCMomentum.cpp
 * @brief Source of the implementation of the vector momentum equation in rotating Rayleigh-Benard convection in a cylindrical annulus
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
#include "Equations/Annulus/Boussinesq/BoussinesqRRBCAnnulusVCMomentum.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "Enums/NonDimensional.hpp"
#include "PhysicalOperators/VelocityAdvection.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   BoussinesqRRBCAnnulusVCMomentum::BoussinesqRRBCAnnulusVCMomentum(SharedEquationParameters spEqParams)
      : IVectorEquation(spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   BoussinesqRRBCAnnulusVCMomentum::~BoussinesqRRBCAnnulusVCMomentum()
   {
   }

   void BoussinesqRRBCAnnulusVCMomentum::setCoupling()
   {
      this->defineCoupling(FieldComponents::Spectral::R, CouplingInformation::PROGNOSTIC, 0, true, false);

      this->defineCoupling(FieldComponents::Spectral::THETA, CouplingInformation::PROGNOSTIC, 0, true, false);

      this->defineCoupling(FieldComponents::Spectral::Z, CouplingInformation::PROGNOSTIC, 0, true, false);
   }

   void BoussinesqRRBCAnnulusVCMomentum::setNLComponents()
   {
      this->addNLComponent(FieldComponents::Spectral::R, 0);

      this->addNLComponent(FieldComponents::Spectral::THETA, 0);

      this->addNLComponent(FieldComponents::Spectral::Z, 0);
   }

   void BoussinesqRRBCAnnulusVCMomentum::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const
   {
      // Get Prandtl number
      MHDFloat Pr = this->eqParams().nd(NonDimensional::PRANDTL);

      /// 
      /// Computation of the advection:
      ///   \f$ \left(\vec u\cdot\nabla\right)u_x\f$
      ///
      FieldComponents::Spectral::Id specId = FieldComponents::Spectral::NOTUSED;
      if(id == FieldComponents::Physical::R)
      {
         specId = FieldComponents::Spectral::R;
      } else if(id == FieldComponents::Physical::THETA)
      {
         specId = FieldComponents::Spectral::THETA;
      } else if(id == FieldComponents::Physical::Z)
      {
         specId = FieldComponents::Spectral::Z;
      }

      Physical::VelocityAdvection<FieldComponents::Physical::R,FieldComponents::Physical::THETA,FieldComponents::Physical::Z>::set(rNLComp, this->unknown().dom(0).phys(), this->unknown().dom(0).grad(specId), 1.0/Pr);
   }

   void BoussinesqRRBCAnnulusVCMomentum::setRequirements()
   {
      // Set temperatur as equation unknown
      this->setName(PhysicalNames::VELOCITY);

      // Set solver timing
      this->setSolveTiming(SolveTiming::PROGNOSTIC);

      // Add velocity to requirements: is scalar?, need spectral?, need physical?, need diff?(, need curl?)
      this->mRequirements.addField(PhysicalNames::VELOCITY, FieldRequirement(false, true, true, true));
   }

}
}
