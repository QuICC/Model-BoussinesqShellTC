/** 
 * @file BoussinesqRRBAnnulusVCMomentum.cpp
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
#include "Equations/Annulus/Boussinesq/BoussinesqRRBAnnulusVCMomentum.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "Enums/NonDimensional.hpp"
#include "PhysicalOperators/VelocityAdvection.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   BoussinesqRRBAnnulusVCMomentum::BoussinesqRRBAnnulusVCMomentum(SharedEquationParameters spEqParams)
      : IVectorEquation(spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   BoussinesqRRBAnnulusVCMomentum::~BoussinesqRRBAnnulusVCMomentum()
   {
   }

   void BoussinesqRRBAnnulusVCMomentum::setCoupling()
   {
      this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::PROGNOSTIC, 0, true, true, false);
   }

   void BoussinesqRRBAnnulusVCMomentum::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const
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

      Physical::VelocityAdvection<FieldComponents::Physical::R,FieldComponents::Physical::THETA,FieldComponents::Physical::Z>::set(rNLComp, this->unknown().dom(0).phys().comp(FieldComponents::Physical::R), this->unknown().dom(0).phys().comp(FieldComponents::Physical::THETA), this->unknown().dom(0).phys().comp(FieldComponents::Physical::Z), this->unknown().dom(0).grad(specId), 1.0/Pr);
   }

   void BoussinesqRRBAnnulusVCMomentum::setRequirements()
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
