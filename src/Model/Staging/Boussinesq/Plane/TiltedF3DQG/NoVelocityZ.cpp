/** 
 * @file NoVelocityZ.cpp
 * @brief Source of the implementation of the non orthogonal vertical velocity computation in the tilted F-plane 3DQG model
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

/// Define small macros allowing to convert to string
#define MAKE_STR_X( _P ) # _P
#define MAKE_STR( _P ) MAKE_STR_X( _P )

// Configuration includes
//

// System includes
//

// External includes

// Class include
//
#include MAKE_STR( QUICC_MODEL_PATH/Boussinesq/Plane/TiltedF3DQG/NoVelocityZ.hpp )

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "Enums/NonDimensional.hpp"
#include "PhysicalOperators/StreamAdvection.hpp"

namespace QuICC {

namespace Equations {

namespace Boussinesq {

namespace Plane {

namespace TiltedF3DQG {

   NoVelocityZ::NoVelocityZ(SharedEquationParameters spEqParams, const SolveTiming::Id time)
      : IScalarEquation(spEqParams)
   {
      // Set solver timing
      this->setSolveTiming(time);

      // Set the variable requirements
      this->setRequirements();
   }

   NoVelocityZ::~NoVelocityZ()
   {
   }

   void NoVelocityZ::setCoupling()
   {
      if(this->solveTiming() == SolveTiming::BEFORE)
      {
         this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::TRIVIAL, 1, true, false, false);
      } else
      {
         this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::TRIVIAL, 1, false, false);
      }
   }

   void NoVelocityZ::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const
   {
      // Assert on scalar component is used
      assert(id == FieldComponents::Physical::SCALAR);

      // Get parameters
      MHDFloat eta3 = std::cos((Math::PI/180.)*this->eqParams().nd(NonDimensional::THETA));

      /// 
      /// Computation of the jacobian:
      ///   \f$ \left(\nabla^{\perp}\psi\cdot\nabla_{\perp}\right)w\f$
      ///
      Physical::StreamAdvection<FieldComponents::Physical::X,FieldComponents::Physical::Y>::set(rNLComp, this->scalar(PhysicalNames::NO_STREAMFUNCTION).dom(0).grad(), this->unknown().dom(0).grad(), 1.0/eta3);
   }

   void NoVelocityZ::setRequirements()
   {
      // Set streamfunction as equation unknown
      this->setName(PhysicalNames::NO_VELOCITYZ);

      // Set non orthogonal vertical velocity requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::NO_VELOCITYZ, FieldRequirement(true, true, true, true));

      // Set non orthogonal streamfunction requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::NO_STREAMFUNCTION, FieldRequirement(true, true, true, true));

      // Gradient does not require Z component
      ArrayB   comps = ArrayB::Constant(3, true);
      comps(0) = false;
      std::map<FieldComponents::Spectral::Id,ArrayB>  gradComps;
      gradComps.insert(std::make_pair(FieldComponents::Spectral::SCALAR, comps));

      // Update temperature gradient requirements
      this->updateFieldRequirements(PhysicalNames::NO_VELOCITYZ).updateGradient(gradComps);

      // Update streamfunction gradient requirements
      this->updateFieldRequirements(PhysicalNames::NO_STREAMFUNCTION).updateGradient(gradComps);
   }

}
}
}
}
}
