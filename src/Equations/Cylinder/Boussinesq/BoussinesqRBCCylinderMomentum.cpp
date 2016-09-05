/** 
 * @file BoussinesqRBCCylinderMomentum.cpp
 * @brief Source of the implementation of the vector momentum equation in Rayleigh-Benard convection in a cylinder (toroidal-poloidal formulation)
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
#include "Equations/Cylinder/Boussinesq/BoussinesqRBCCylinderMomentum.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "Enums/NonDimensional.hpp"
#include "PhysicalOperators/VelocityAdvection.hpp"
#include "PhysicalOperators/Cross.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   BoussinesqRBCCylinderMomentum::BoussinesqRBCCylinderMomentum(SharedEquationParameters spEqParams)
      : IVectorEquation(spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   BoussinesqRBCCylinderMomentum::~BoussinesqRBCCylinderMomentum()
   {
   }

   void BoussinesqRBCCylinderMomentum::setCoupling()
   {
      this->defineCoupling(FieldComponents::Spectral::TOR, CouplingInformation::PROGNOSTIC, 0, true, false);

      this->defineCoupling(FieldComponents::Spectral::POL, CouplingInformation::PROGNOSTIC, 0, true, false);
   }

   void BoussinesqRBCCylinderMomentum::setNLComponents()
   {
      this->addNLComponent(FieldComponents::Spectral::TOR, 0);

      this->addNLComponent(FieldComponents::Spectral::POL, 0);
   }

   void BoussinesqRBCCylinderMomentum::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id compId) const
   {
      // Get Prandtl number
      MHDFloat Pr = this->eqParams().nd(NonDimensional::PRANDTL);

      ///
      /// Compute \f$\left(\nabla\wedge\vec u\right)\wedge\vec u\f$
      ///
      switch(compId)
      {
         case(FieldComponents::Physical::R):
            Physical::Cross<FieldComponents::Physical::THETA,FieldComponents::Physical::Z>::set(rNLComp, this->unknown().dom(0).curl(), this->unknown().dom(0).phys(), 1.0/Pr);
            break;
         case(FieldComponents::Physical::THETA):
            Physical::Cross<FieldComponents::Physical::Z,FieldComponents::Physical::R>::set(rNLComp, this->unknown().dom(0).curl(), this->unknown().dom(0).phys(), 1.0/Pr);
            break;
         case(FieldComponents::Physical::Z):
            Physical::Cross<FieldComponents::Physical::R,FieldComponents::Physical::THETA>::set(rNLComp, this->unknown().dom(0).curl(), this->unknown().dom(0).phys(), 1.0/Pr);
            break;
         default:
            assert(false);
            break;
      }
   }

   void BoussinesqRBCCylinderMomentum::setRequirements()
   {
      // Set temperatur as equation unknown
      this->setName(PhysicalNames::VELOCITY);

      // Set solver timing
      this->setSolveTiming(SolveTiming::PROGNOSTIC);

      // Add X velocity to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::VELOCITY, FieldRequirement(false, true, true, true));
   }

}
}
