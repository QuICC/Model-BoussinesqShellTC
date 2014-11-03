/** 
 * @file BoussinesqFPlane3DQGVelocityZ.cpp
 * @brief Source of the implementation of the upright vertical velocity equation in the F-plane 3DQG model
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
#include "Equations/Asymptotics/FPlane3DQG/Boussinesq/BoussinesqFPlane3DQGVelocityZ.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "PhysicalOperators/StreamAdvection.hpp"
#include "TypeSelectors/EquationEigenSelector.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   BoussinesqFPlane3DQGVelocityZ::BoussinesqFPlane3DQGVelocityZ(SharedEquationParameters spEqParams)
      : IScalarEquation(spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   BoussinesqFPlane3DQGVelocityZ::~BoussinesqFPlane3DQGVelocityZ()
   {
   }

   void BoussinesqFPlane3DQGVelocityZ::setCoupling()
   {
      this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::PROGNOSTIC, 1, true, true, false);
   }

   void BoussinesqFPlane3DQGVelocityZ::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const
   {
      // Assert on scalar component is used
      assert(id == FieldComponents::Physical::SCALAR);

      /// 
      /// Computation of the jacobian:
      ///   \f$ \left(\nabla^{\perp}\psi\cdot\nabla_{\perp}\right)w\f$
      ///
      Physical::StreamAdvection<FieldComponents::Physical::TWO,FieldComponents::Physical::THREE>::set(rNLComp, this->scalar(PhysicalNames::STREAMFUNCTION).dom(0).grad(), this->unknown().dom(0).grad(), 1.0);
   }

   void BoussinesqFPlane3DQGVelocityZ::setRequirements()
   {
      // Set vertical velocity as equation unknown
      this->setName(PhysicalNames::VELOCITYZ);

      // Set solver timing
      this->setSolveTiming(SolveTiming::PROGNOSTIC);

      // Add vertical velocity requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::VELOCITYZ, FieldRequirement(true, true, true, true));

      // Set non orthogonal streamfunction requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::STREAMFUNCTION, FieldRequirement(true, true, true, true));
   }

}
}
