/** 
 * @file BoussinesqPrecessionSphereMomentum.cpp
 * @brief Source of the implementation of the vector Navier-Stokes equation in the Boussinesq precession in a sphere model
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Configuration includes
//
#include "TypeSelectors/TransformSelector.hpp"

// System includes
//

// External includes
//

// Class include
//
#include "Equations/Sphere/Boussinesq/BoussinesqPrecessionSphereMomentum.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "Enums/NonDimensional.hpp"
#include "PhysicalOperators/Cross.hpp"
#include "PhysicalOperators/SphericalPrecession.hpp"
#include "PhysicalOperators/SphericalPoincare.hpp"

namespace QuICC {

namespace Equations {

   BoussinesqPrecessionSphereMomentum::BoussinesqPrecessionSphereMomentum(SharedEquationParameters spEqParams)
      : IVectorEquation(spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   BoussinesqPrecessionSphereMomentum::~BoussinesqPrecessionSphereMomentum()
   {
   }

   void BoussinesqPrecessionSphereMomentum::setCoupling()
   {
      #if defined QUICC_SPATIALSCHEME_BLFL || defined QUICC_SPATIALSCHEME_WLFL
         int start = 1;
      #elif defined QUICC_SPATIALSCHEME_BLFM || defined QUICC_SPATIALSCHEME_WLFM
         int start = 0;
      #endif //defined QUICC_SPATIALSCHEME_BLFL || defined QUICC_SPATIALSCHEME_WLFL

      this->defineCoupling(FieldComponents::Spectral::TOR, CouplingInformation::PROGNOSTIC, start, true, false);

      this->defineCoupling(FieldComponents::Spectral::POL, CouplingInformation::PROGNOSTIC, start, true, false);

      // Create R, theta and phi physical grids
      int nR = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D,Dimensions::Space::PHYSICAL);
      this->mR = Transform::TransformSelector<Dimensions::Transform::TRA1D>::Type::generateGrid(nR);
      int nTh = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM2D,Dimensions::Space::PHYSICAL);
      this->mTheta = Transform::TransformSelector<Dimensions::Transform::TRA2D>::Type::generateGrid(nTh);
      int nPh = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM3D,Dimensions::Space::PHYSICAL);
      this->mPhi = Transform::TransformSelector<Dimensions::Transform::TRA3D>::Type::generateGrid(nPh);
   }

   void BoussinesqPrecessionSphereMomentum::setNLComponents()
   {
      this->addNLComponent(FieldComponents::Spectral::TOR, 0);

      this->addNLComponent(FieldComponents::Spectral::POL, 0);
   }

   void BoussinesqPrecessionSphereMomentum::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id compId) const
   {  
      ///
      /// Compute \f$\left(\nabla\wedge\vec u\right)\wedge\vec u\f$
      ///
      switch(compId)
      {
         case(FieldComponents::Physical::R):
            Physical::Cross<FieldComponents::Physical::THETA,FieldComponents::Physical::PHI>::set(rNLComp, this->unknown().dom(0).curl(), this->unknown().dom(0).phys(), 1.0);
            break;
         case(FieldComponents::Physical::THETA):
            Physical::Cross<FieldComponents::Physical::PHI,FieldComponents::Physical::R>::set(rNLComp, this->unknown().dom(0).curl(), this->unknown().dom(0).phys(), 1.0);
            break;
         case(FieldComponents::Physical::PHI):
            Physical::Cross<FieldComponents::Physical::R,FieldComponents::Physical::THETA>::set(rNLComp, this->unknown().dom(0).curl(), this->unknown().dom(0).phys(), 1.0);
            break;
         default:
            assert(false);
            break;
      }

      ///
      /// Compute Coriolis + Precession term
      ///
      MHDFloat Po = this->eqParams().nd(NonDimensional::POINCARE);
      MHDFloat alpha = (this->eqParams().nd(NonDimensional::ALPHA)*Math::PI/180.);
      #if defined QUICC_SPATIALSCHEME_BLFL || defined QUICC_SPATIALSCHEME_WLFL
         MHDFloat corC = 1.0;
      #elif defined QUICC_SPATIALSCHEME_BLFM || defined QUICC_SPATIALSCHEME_WLFM
         MHDFloat corC = 0.0;
      #endif //defined QUICC_SPATIALSCHEME_BLFL || defined QUICC_SPATIALSCHEME_WLFL
      Physical::SphericalPrecession::add(rNLComp, compId, this->unknown().dom(0).spRes(), this->mTheta, this->mPhi, this->unknown().dom(0).phys(), this->time(), alpha, corC, Po, 2.0);

      /// Compute Poincare term
      Physical::SphericalPoincare::add(rNLComp, compId, this->unknown().dom(0).spRes(), this->mR, this->mTheta, this->mPhi, this->unknown().dom(0).phys(), this->time(), alpha, Po);
   }

   void BoussinesqPrecessionSphereMomentum::setRequirements()
   {
      // Set velocity as equation unknown
      this->setName(PhysicalNames::VELOCITY);

      // Set solver timing
      this->setSolveTiming(SolveTiming::PROGNOSTIC);

      // Add velocity to requirements: is scalar?, need spectral?, need physical?, need diff?(, need curl?)
      this->mRequirements.addField(PhysicalNames::VELOCITY, FieldRequirement(false, true, true, false, true));
   }

}
}
