/** 
 * @file BoussinesqRBCAnnulusVCTransport.cpp
 * @brief Source of the implementation of the transport equation in Rayleigh-Benard convection in a cylindrical annulus
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
#include "Equations/Annulus/Boussinesq/BoussinesqRBCAnnulusVCTransport.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "Enums/NonDimensional.hpp"
#include "PhysicalOperators/VelocityHeatAdvection.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   BoussinesqRBCAnnulusVCTransport::BoussinesqRBCAnnulusVCTransport(SharedEquationParameters spEqParams)
      : IScalarEquation(spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   BoussinesqRBCAnnulusVCTransport::~BoussinesqRBCAnnulusVCTransport()
   {
   }

   void BoussinesqRBCAnnulusVCTransport::setCoupling()
   {
      this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::PROGNOSTIC, 0, true, false);
   }

   void BoussinesqRBCAnnulusVCTransport::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const
   {
      // Assert on scalar component is used
      assert(id == FieldComponents::Physical::SCALAR);

      /// 
      /// Computation of the advection:
      ///   \f$ \left(\vec u\cdot\nabla\right)\theta\f$
      ///
      Physical::VelocityHeatAdvection<FieldComponents::Physical::R,FieldComponents::Physical::THETA,FieldComponents::Physical::Z>::set(rNLComp, this->vector(PhysicalNames::VELOCITY).dom(0).phys(), this->unknown().dom(0).grad(), 1.0);
   }

   void BoussinesqRBCAnnulusVCTransport::setRequirements()
   {
      // Set temperatur as equation unknown
      this->setName(PhysicalNames::TEMPERATURE);

      // Set solver timing
      this->setSolveTiming(SolveTiming::PROGNOSTIC);

      // Add temperature to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::TEMPERATURE, FieldRequirement(true, true, false, true));

      // Add velocity to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::VELOCITY, FieldRequirement(false, true, true, false));
   }

}
}
