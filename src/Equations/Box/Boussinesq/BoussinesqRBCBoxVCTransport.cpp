/** 
 * @file BoussinesqRBCBoxVCTransport.cpp
 * @brief Source of the implementation of the transport equation in Rayleigh-Benard convection in 3D box
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
#include "Equations/Box/Boussinesq/BoussinesqRBCBoxVCTransport.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "Enums/NonDimensional.hpp"
#include "PhysicalOperators/VelocityHeatAdvection.hpp"

namespace QuICC {

namespace Equations {

   BoussinesqRBCBoxVCTransport::BoussinesqRBCBoxVCTransport(SharedEquationParameters spEqParams)
      : IScalarEquation(spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   BoussinesqRBCBoxVCTransport::~BoussinesqRBCBoxVCTransport()
   {
   }

   void BoussinesqRBCBoxVCTransport::setCoupling()
   {
      this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::PROGNOSTIC, 0, true, false);
   }

   void BoussinesqRBCBoxVCTransport::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const
   {
      // Assert on scalar component is used
      assert(id == FieldComponents::Physical::SCALAR);

      /// 
      /// Computation of the advection:
      ///   \f$ \left(\vec u\cdot\nabla\right)\theta\f$
      ///
      Physical::VelocityHeatAdvection<FieldComponents::Physical::X,FieldComponents::Physical::Y,FieldComponents::Physical::Z>::set(rNLComp, this->vector(PhysicalNames::VELOCITY).dom(0).phys(), this->unknown().dom(0).grad(), 1.0);
   }

   void BoussinesqRBCBoxVCTransport::setRequirements()
   {
      // Set temperatur as equation unknown
      this->setName(PhysicalNames::TEMPERATURE);

      // Set solver timing
      this->setSolveTiming(SolveTiming::PROGNOSTIC);

      // Add X velocity to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::TEMPERATURE, FieldRequirement(true, true, false, true));

      // Add X velocity to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::VELOCITY, FieldRequirement(false, true, true, false));
   }

}
}
