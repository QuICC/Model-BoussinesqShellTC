/** 
 * @file BoussinesqRBCPlane2DVCMomentum.cpp
 * @brief Source of the implementation of the vector momentum equation in Rayleigh-Benard convection in a plane layer (2D)
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
#include "Equations/Box/Boussinesq/BoussinesqRBCPlane2DVCMomentum.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "Enums/NonDimensional.hpp"
#include "PhysicalOperators/VelocityAdvection.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   BoussinesqRBCPlane2DVCMomentum::BoussinesqRBCPlane2DVCMomentum(SharedEquationParameters spEqParams)
      : IVectorEquation(spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   BoussinesqRBCPlane2DVCMomentum::~BoussinesqRBCPlane2DVCMomentum()
   {
   }

   void BoussinesqRBCPlane2DVCMomentum::setCoupling()
   {
      this->defineCoupling(FieldComponents::Spectral::X, CouplingInformation::PROGNOSTIC, 1, true, false);

      this->defineCoupling(FieldComponents::Spectral::Z, CouplingInformation::PROGNOSTIC, 1, true, false);
   }

   void BoussinesqRBCPlane2DVCMomentum::setNLComponents()
   {
      this->addNLComponent(FieldComponents::Spectral::X, 0);

      this->addNLComponent(FieldComponents::Spectral::Z, 0);
   }

   void BoussinesqRBCPlane2DVCMomentum::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const
   {
      // Get Prandtl number
      MHDFloat Pr = this->eqParams().nd(NonDimensional::PRANDTL);

      /// 
      /// Computation of the advection:
      ///   \f$ \left(\vec u\cdot\nabla\right)u_x\f$
      ///
      FieldComponents::Spectral::Id specId = FieldComponents::Spectral::NOTUSED;
      if(id == FieldComponents::Physical::X)
      {
         specId = FieldComponents::Spectral::X;
      } else if(id == FieldComponents::Physical::Z)
      {
         specId = FieldComponents::Spectral::Z;
      }

      Physical::VelocityAdvection<FieldComponents::Physical::X,FieldComponents::Physical::Z>::set(rNLComp, this->unknown().dom(0).phys(), this->unknown().dom(0).grad(specId), 1.0/Pr);
   }

   void BoussinesqRBCPlane2DVCMomentum::setRequirements()
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
