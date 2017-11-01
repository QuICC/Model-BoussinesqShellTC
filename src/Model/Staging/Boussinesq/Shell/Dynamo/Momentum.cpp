/** 
 * @file Momentum.cpp
 * @brief Source of the implementation of the vector Navier-Stokes equation in the Boussinesq thermal convection dynamo in a spherical shell model
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

/// Define small macros allowing to convert to string
#define MAKE_STR_X( _P ) # _P
#define MAKE_STR( _P ) MAKE_STR_X( _P )

// Configuration includes
//
#include "TypeSelectors/TransformSelector.hpp"

// System includes
//

// External includes
//

// Class include
//
#include MAKE_STR( QUICC_MODEL_PATH/Boussinesq/Shell/Dynamo/Momentum.hpp )

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "Enums/NonDimensional.hpp"
#include "PhysicalOperators/Cross.hpp"
#include "PhysicalOperators/SphericalCoriolis.hpp"

namespace QuICC {

namespace Equations {

namespace Boussinesq {

namespace Shell {

namespace Dynamo {

   Momentum::Momentum(SharedEquationParameters spEqParams)
      : IVectorEquation(spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   Momentum::~Momentum()
   {
   }

   void Momentum::setCoupling()
   {
      #ifdef QUICC_SPATIALSCHEME_SLFL
         int start = 1;
      #else //if QUICC_SPATIALSCHEME_SLFM
         int start = 0;
      #endif //QUICC_SPATIALSCHEME_SLFL

      this->defineCoupling(FieldComponents::Spectral::TOR, CouplingInformation::PROGNOSTIC, start, true, false);

      this->defineCoupling(FieldComponents::Spectral::POL, CouplingInformation::PROGNOSTIC, start, true, false);

      #ifdef QUICC_SPATIALSCHEME_SLFL
         // Create cos(theta) and sin(theta) data for Coriolis term
         int nTh = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM2D,Dimensions::Space::PHYSICAL);
         Array thGrid = Transform::TransformSelector<Dimensions::Transform::TRA2D>::Type::generateGrid(nTh);
         this->mCosTheta = thGrid.array().cos();
         this->mSinTheta = thGrid.array().sin();
      #endif //QUICC_SPATIALSCHEME_SLFL
   }

   void Momentum::setNLComponents()
   {
      this->addNLComponent(FieldComponents::Spectral::TOR, 0);

      this->addNLComponent(FieldComponents::Spectral::POL, 0);
   }

   void Momentum::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id compId) const
   {
      // Get square root of Taylor number
      MHDFloat T = std::sqrt(this->eqParams().nd(NonDimensional::TAYLOR));
      MHDFloat Pm = this->eqParams().nd(NonDimensional::MAGPRANDTL);

      ///
      /// Compute \f$\vec u\wedge\left(\nabla\wedge\vec u\right) + \left(\nabla\wedge\vec B\right)\wedge\vec B\f$
      ///
      switch(compId)
      {
         case(FieldComponents::Physical::R):
            Physical::Cross<FieldComponents::Physical::THETA,FieldComponents::Physical::PHI>::set(rNLComp, this->unknown().dom(0).curl(), this->unknown().dom(0).phys(), 1.0);
            Physical::Cross<FieldComponents::Physical::THETA,FieldComponents::Physical::PHI>::add(rNLComp, this->vector(PhysicalNames::MAGNETIC).dom(0).phys(), this->vector(PhysicalNames::MAGNETIC).dom(0).curl(), T/Pm);
            break;
         case(FieldComponents::Physical::THETA):
            Physical::Cross<FieldComponents::Physical::PHI,FieldComponents::Physical::R>::set(rNLComp, this->unknown().dom(0).curl(), this->unknown().dom(0).phys(), 1.0);
            Physical::Cross<FieldComponents::Physical::PHI,FieldComponents::Physical::R>::add(rNLComp, this->vector(PhysicalNames::MAGNETIC).dom(0).phys(), this->vector(PhysicalNames::MAGNETIC).dom(0).curl(), T/Pm);
            break;
         case(FieldComponents::Physical::PHI):
            Physical::Cross<FieldComponents::Physical::R,FieldComponents::Physical::THETA>::set(rNLComp, this->unknown().dom(0).curl(), this->unknown().dom(0).phys(), 1.0);
            Physical::Cross<FieldComponents::Physical::R,FieldComponents::Physical::THETA>::add(rNLComp, this->vector(PhysicalNames::MAGNETIC).dom(0).phys(), this->vector(PhysicalNames::MAGNETIC).dom(0).curl(), T/Pm);
            break;
         default:
            assert(false);
            break;
      }

      #ifdef QUICC_SPATIALSCHEME_SLFL
         ///
         /// Compute Coriolis term
         ///
         Physical::SphericalCoriolis::add(rNLComp, compId, this->unknown().dom(0).spRes(), this->mCosTheta, this->mSinTheta, this->unknown().dom(0).phys(), T);
      #endif //QUICC_SPATIALSCHEME_SLFL
   }

   void Momentum::setRequirements()
   {
      // Set velocity as equation unknown
      this->setName(PhysicalNames::VELOCITY);

      // Set solver timing
      this->setSolveTiming(SolveTiming::PROGNOSTIC);

      // Add velocity to requirements: is scalar?, need spectral?, need physical?, need diff?(, need curl?)
      this->mRequirements.addField(PhysicalNames::VELOCITY, FieldRequirement(false, true, true, false, true));

      // Add magnetic to requirements: is scalar?, need spectral?, need physical?, need diff?(, need curl?)
      this->mRequirements.addField(PhysicalNames::MAGNETIC, FieldRequirement(false, true, true, false, true));
   }

}
}
}
}
}