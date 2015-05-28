/** 
 * @file BoussinesqBeta3DQGVelocityZ.cpp
 * @brief Source of the implementation of the vertical velocity computation in the Beta 3DQG model
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Configuration includes
//

// System includes
//

// External includes

// Class include
//
#include "Equations/Asymptotics/Beta3DQG/Boussinesq/BoussinesqBeta3DQGVelocityZ.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "Enums/NonDimensional.hpp"
#include "PhysicalOperators/StreamAdvection.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   BoussinesqBeta3DQGVelocityZ::BoussinesqBeta3DQGVelocityZ(SharedEquationParameters spEqParams)
      : IScalarEquation(spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   BoussinesqBeta3DQGVelocityZ::~BoussinesqBeta3DQGVelocityZ()
   {
   }

   void BoussinesqBeta3DQGVelocityZ::setCoupling()
   {
      this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::PROGNOSTIC, 0, true, false);
   }

   void BoussinesqBeta3DQGVelocityZ::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const
   {
      // Assert on scalar component is used
      assert(id == FieldComponents::Physical::SCALAR);

      /// 
      /// Computation of the jacobian:
      ///   \f$ \left(\nabla^{\perp}\psi\cdot\nabla_{\perp}\right)w\f$
      ///
      Physical::StreamAdvection<>::set(rNLComp, this->scalar(PhysicalNames::STREAMFUNCTION).dom(0).grad(), this->unknown().dom(0).grad(), 0.0);
   }

   void BoussinesqBeta3DQGVelocityZ::setRequirements()
   {
      // Set streamfunction as equation unknown
      this->setName(PhysicalNames::VELOCITYZ);

      // Set solver timing
      this->setSolveTiming(SolveTiming::PROGNOSTIC);

      // Set vorticity requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::VELOCITYZ, FieldRequirement(true, true, false, true));

      // Add streamfunction requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::STREAMFUNCTION, FieldRequirement(true, true, false, true));
   }

}
}
