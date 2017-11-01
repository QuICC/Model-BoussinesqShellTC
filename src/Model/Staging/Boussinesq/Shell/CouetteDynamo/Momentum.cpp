/** 
 * @file Momentum.cpp
 * @brief Source of the implementation of the vector Navier-Stokes equation in the Boussinesq spherical Couette dynamo in a spherical shell model
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
#include MAKE_STR( QUICC_MODEL_PATH/Boussinesq/Shell/CouetteDynamo/Momentum.hpp )

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

namespace CouetteDynamo {

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

      this->defineCoupling(FieldComponents::Spectral::TOR, CouplingInformation::PROGNOSTIC, start, true, false, true);

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
      MHDFloat Ro = std::abs(this->eqParams().nd(NonDimensional::ROSSBY));
      MHDFloat Lambda = this->eqParams().nd(NonDimensional::MODELSASSER);
      ///
      /// Compute \f$\left(\nabla\wedge\vec u\right)\wedge\vec u\f$
      ///
      switch(compId)
      {
      // TODO: ask Philippe if there is a way to compute \hat{z}+b on spectral space and therefore be more efficient
         case(FieldComponents::Physical::R):
            Physical::Cross<FieldComponents::Physical::THETA,FieldComponents::Physical::PHI>::set(rNLComp, this->unknown().dom(0).curl(), this->unknown().dom(0).phys(), Ro);
            Physical::Cross<FieldComponents::Physical::THETA,FieldComponents::Physical::PHI>::add(rNLComp, this->vector(PhysicalNames::MAGNETIC).dom(0).phys(), this->vector(PhysicalNames::MAGNETIC).dom(0).curl(), Lambda);
            Physical::Cross<FieldComponents::Physical::THETA,FieldComponents::Physical::PHI>::add(rNLComp, this->vector(PhysicalNames::IMPOSED_MAGNETIC).dom(0).phys(), this->vector(PhysicalNames::MAGNETIC).dom(0).curl(), Lambda);
            break;
         case(FieldComponents::Physical::THETA):
            Physical::Cross<FieldComponents::Physical::PHI,FieldComponents::Physical::R>::set(rNLComp, this->unknown().dom(0).curl(), this->unknown().dom(0).phys(), Ro);
            Physical::Cross<FieldComponents::Physical::PHI,FieldComponents::Physical::R>::add(rNLComp, this->vector(PhysicalNames::MAGNETIC).dom(0).phys(), this->vector(PhysicalNames::MAGNETIC).dom(0).curl(), Lambda);
            Physical::Cross<FieldComponents::Physical::PHI,FieldComponents::Physical::R>::add(rNLComp, this->vector(PhysicalNames::IMPOSED_MAGNETIC).dom(0).phys(), this->vector(PhysicalNames::MAGNETIC).dom(0).curl(), Lambda);
            break;
         case(FieldComponents::Physical::PHI):
            Physical::Cross<FieldComponents::Physical::R,FieldComponents::Physical::THETA>::set(rNLComp, this->unknown().dom(0).curl(), this->unknown().dom(0).phys(), Ro);
            Physical::Cross<FieldComponents::Physical::R,FieldComponents::Physical::THETA>::add(rNLComp, this->vector(PhysicalNames::MAGNETIC).dom(0).phys(), this->vector(PhysicalNames::MAGNETIC).dom(0).curl(), Lambda);
            Physical::Cross<FieldComponents::Physical::R,FieldComponents::Physical::THETA>::add(rNLComp, this->vector(PhysicalNames::IMPOSED_MAGNETIC).dom(0).phys(), this->vector(PhysicalNames::MAGNETIC).dom(0).curl(), Lambda);
            break;
         default:
            assert(false);
            break;
      }

      #ifdef QUICC_SPATIALSCHEME_SLFL
         ///
         /// Compute Coriolis term
         ///
         Physical::SphericalCoriolis::add(rNLComp, compId, this->unknown().dom(0).spRes(), this->mCosTheta, this->mSinTheta, this->unknown().dom(0).phys(), 2.0);
      #endif //QUICC_SPATIALSCHEME_SLFL
   }

   Datatypes::SpectralScalarType::PointType Momentum::boundaryValue(FieldComponents::Spectral::Id compId, const int i, const int j, const int k) const
   {  
      if(compId == FieldComponents::Spectral::TOR)
      {
         #if defined QUICC_SPATIALSCHEME_SLFL
            if(this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k) == 1)
            {
               if(this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j,k) == 0)
               {
                  if(i == 1)
                  {
         #elif defined QUICC_SPATIALSCHEME_SLFM
            if(this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k) == 0)
            {
               if(this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j,k) == 1)
               {
                  if(i == 1)
                  {
         #endif //defined QUICC_SPATIALSCHEME_SLFL

                     MHDFloat sgn = std::pow(-1.0,std::signbit(this->eqParams().nd(NonDimensional::ROSSBY)));
                     MHDFloat ri = this->eqParams().nd(NonDimensional::RO)*this->eqParams().nd(NonDimensional::RRATIO);
                     MHDFloat norm = std::sqrt(3.0/(4.0*Math::PI));
                     return Datatypes::SpectralScalarType::PointType(sgn*ri/norm);
                  }
               }
            }

         return Datatypes::SpectralScalarType::PointType(0.0);
      } else
      {
         throw std::logic_error("Poloidal component should not set inhomogeneous boundary condition");

         return Datatypes::SpectralScalarType::PointType(0.0);
      }
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

      // Add imposed magnetic field to requirements: is scalar?, need spectral?, need physical?, need diff?(, need curl?)
      this->mImposedRequirements.addField(PhysicalNames::IMPOSED_MAGNETIC, FieldRequirement(false, true, true, false, false));
   }

}
}
}
}
}
