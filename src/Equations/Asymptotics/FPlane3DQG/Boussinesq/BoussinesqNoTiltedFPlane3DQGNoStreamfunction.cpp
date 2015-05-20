/** 
 * @file BoussinesqNoTiltedFPlane3DQGNoStreamfunction.cpp
 * @brief Source of the implementation of the non orthogonal streamfuntion computation in the tilted F-plane 3DQG model
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Configuration includes
//

// System includes
//

// External includes

// Class include
//
#include "Equations/Asymptotics/FPlane3DQG/Boussinesq/BoussinesqNoTiltedFPlane3DQGNoStreamfunction.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "Enums/NonDimensional.hpp"
#include "PhysicalOperators/StreamAdvection.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   BoussinesqNoTiltedFPlane3DQGNoStreamfunction::BoussinesqNoTiltedFPlane3DQGNoStreamfunction(SharedEquationParameters spEqParams)
      : IScalarEquation(spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   BoussinesqNoTiltedFPlane3DQGNoStreamfunction::~BoussinesqNoTiltedFPlane3DQGNoStreamfunction()
   {
   }

   void BoussinesqNoTiltedFPlane3DQGNoStreamfunction::setCoupling()
   {
      this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::PROGNOSTIC, 1, true, false);
   }

   void BoussinesqNoTiltedFPlane3DQGNoStreamfunction::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const
   {
      // Assert on scalar component is used
      assert(id == FieldComponents::Physical::SCALAR);

      // Get paramters
      MHDFloat eta3 = std::cos((Math::PI/180.)*this->eqParams().nd(NonDimensional::THETA));

      /// 
      /// Computation of the jacobian:
      ///   \f$ \left(\nabla^{\perp}\psi\cdot\nabla_{\perp}\right)\nabla^2_{\perp}\psi\f$
      ///
      Physical::StreamAdvection<FieldComponents::Physical::TWO,FieldComponents::Physical::THREE>::set(rNLComp, this->unknown().dom(0).grad(), this->scalar(PhysicalNames::NO_VORTICITYZ).dom(0).grad(), 1.0/eta3);
   }

   void BoussinesqNoTiltedFPlane3DQGNoStreamfunction::setRequirements()
   {
      // Set streamfunction as equation unknown
      this->setName(PhysicalNames::NO_STREAMFUNCTION);

      // Set solver timing
      this->setSolveTiming(SolveTiming::PROGNOSTIC);

      // Set non orthogonal streamfunction requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::NO_STREAMFUNCTION, FieldRequirement(true, true, true, true));

      // Set non orthogonal vertical vorticity requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::NO_VORTICITYZ, FieldRequirement(true, true, true, true));
   }

}
}
