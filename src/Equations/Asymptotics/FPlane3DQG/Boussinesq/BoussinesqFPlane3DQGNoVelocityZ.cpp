/** 
 * @file BoussinesqFPlane3DQGNoVelocityZ.cpp
 * @brief Source of the implementation of the non orthogonal vertical velocity computation in the F-plane 3DQG model
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Configuration includes
//

// System includes
//

// External includes

// Class include
//
#include "Equations/Asymptotics/FPlane3DQG/Boussinesq/BoussinesqFPlane3DQGNoVelocityZ.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "Enums/NonDimensional.hpp"
#include "PhysicalOperators/StreamAdvection.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   BoussinesqFPlane3DQGNoVelocityZ::BoussinesqFPlane3DQGNoVelocityZ(SharedEquationParameters spEqParams, const SolveTiming::Id time)
      : IScalarEquation(spEqParams)
   {
      // Set solver timing
      this->setSolveTiming(time);

      // Set the variable requirements
      this->setRequirements();
   }

   BoussinesqFPlane3DQGNoVelocityZ::~BoussinesqFPlane3DQGNoVelocityZ()
   {
   }

   void BoussinesqFPlane3DQGNoVelocityZ::setCoupling()
   {
      if(this->solveTiming() == SolveTiming::BEFORE)
      {
         this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::TRIVIAL, 1, true, false, false, false);
      } else
      {
         this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::TRIVIAL, 1, false, false, false);
      }
   }

   void BoussinesqFPlane3DQGNoVelocityZ::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const
   {
      // Assert on scalar component is used
      assert(id == FieldComponents::Physical::SCALAR);

      // Get parameters
      MHDFloat eta3 = std::cos((Math::PI/180.)*this->eqParams().nd(NonDimensional::THETA));

      /// 
      /// Computation of the jacobian:
      ///   \f$ \left(\nabla^{\perp}\psi\cdot\nabla_{\perp}\right)w\f$
      ///
      Physical::StreamAdvection<FieldComponents::Physical::TWO,FieldComponents::Physical::THREE>::set(rNLComp, this->scalar(PhysicalNames::NO_STREAMFUNCTION).dom(0).grad(), this->unknown().dom(0).grad(), 1.0/eta3);
   }

   void BoussinesqFPlane3DQGNoVelocityZ::setRequirements()
   {
      // Set streamfunction as equation unknown
      this->setName(PhysicalNames::NO_VELOCITYZ);

      // Set non orthogonal vertical velocity requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::NO_VELOCITYZ, FieldRequirement(true, true, true, true));

      // Set non orthogonal streamfunction requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::NO_STREAMFUNCTION, FieldRequirement(true, true, true, true));
   }

}
}
