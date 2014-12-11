/** 
 * @file BoussinesqRRB1DBoxVCMomentum.cpp
 * @brief Source of the implementation of the vector momentum equation for rotating Rayleigh-Benard convection in 1D box
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
#include "Equations/Box/Boussinesq/BoussinesqRRB1DBoxVCMomentum.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "Enums/NonDimensional.hpp"
#include "PhysicalOperators/VelocityAdvection.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   BoussinesqRRB1DBoxVCMomentum::BoussinesqRRB1DBoxVCMomentum(SharedEquationParameters spEqParams)
      : IVectorEquation(spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   BoussinesqRRB1DBoxVCMomentum::~BoussinesqRRB1DBoxVCMomentum()
   {
   }

   void BoussinesqRRB1DBoxVCMomentum::setCoupling()
   {
      this->defineCoupling(FieldComponents::Spectral::X, CouplingInformation::PROGNOSTIC, 1, true, true, false);

      this->defineCoupling(FieldComponents::Spectral::Y, CouplingInformation::PROGNOSTIC, 1, true, true, false);

      this->defineCoupling(FieldComponents::Spectral::Z, CouplingInformation::PROGNOSTIC, 1, true, true, false);
   }

   void BoussinesqRRB1DBoxVCMomentum::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const
   {
      /// 
      /// Computation of the advection:
      ///   \f$ \left(\vec u\cdot\nabla\right)u_x\f$
      ///
      FieldComponents::Spectral::Id specId;
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
      
      Physical::VelocityAdvection<FieldComponents::Physical::X,FieldComponents::Physical::Y,FieldComponents::Physical::Z>::set(rNLComp, this->unknown().dom(0).phys().comp(FieldComponents::Physical::X), this->unknown().dom(0).phys().comp(FieldComponents::Physical::Y), this->unknown().dom(0).phys().comp(FieldComponents::Physical::Z), this->unknown().dom(0).grad(specId), 1.0);
   }

   void BoussinesqRRB1DBoxVCMomentum::setRequirements()
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
