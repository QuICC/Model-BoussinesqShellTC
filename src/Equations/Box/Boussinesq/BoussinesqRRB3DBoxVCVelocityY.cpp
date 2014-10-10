/** 
 * @file BoussinesqRRB3DBoxVCVelocityY.cpp
 * @brief Source of the implementation of th momentum equation for the Y component in rotating Rayleigh-Benard convection in 3D box
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
#include "Equations/Box/Boussinesq/BoussinesqRRB3DBoxVCVelocityY.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "Enums/NonDimensional.hpp"
#include "PhysicalOperators/VelocityAdvection.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   BoussinesqRRB3DBoxVCVelocityY::BoussinesqRRB3DBoxVCVelocityY(SharedEquationParameters spEqParams)
      : IScalarEquation(spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   BoussinesqRRB3DBoxVCVelocityY::~BoussinesqRRB3DBoxVCVelocityY()
   {
   }

   void BoussinesqRRB3DBoxVCVelocityY::setCoupling()
   {
      this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::PROGNOSTIC, 0, true, true, false);
   }

   void BoussinesqRRB3DBoxVCVelocityY::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const
   {
      // Assert on scalar component is used
      assert(id == FieldComponents::Physical::SCALAR);

      /// 
      /// Computation of the advection:
      ///   \f$ \left(\vec u\cdot\nabla\right)u_y\f$
      ///
      Physical::VelocityAdvection<FieldComponents::Physical::ONE,FieldComponents::Physical::TWO,FieldComponents::Physical::THREE>::set(rNLComp, this->scalar(PhysicalNames::VELOCITYX).dom(0).phys(), this->scalar(PhysicalNames::VELOCITYY).dom(0).phys(), this->scalar(PhysicalNames::VELOCITYZ).dom(0).phys(), this->unknown().dom(0).grad(), 1.0);
   }

   void BoussinesqRRB3DBoxVCVelocityY::setRequirements()
   {
      // Set temperatur as equation unknown
      this->setName(PhysicalNames::VELOCITYY);

      // Set solver timing
      this->setSolveTiming(SolveTiming::PROGNOSTIC);

      // Add Y velocity to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::VELOCITYY, FieldRequirement(true, true, true, true));

      // Add X velocity to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::VELOCITYX, FieldRequirement(true, true, true, false));

      // Add Z velocity to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::VELOCITYZ, FieldRequirement(true, true, true, false));
   }

}
}
