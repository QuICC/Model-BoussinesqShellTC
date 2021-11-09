/** 
 * @file Momentum.cpp
 * @brief Source of the implementation of the vector Navier-Stokes equation in the Boussinesq rotating thermal convection in a sphere model
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
#include MAKE_STR( QUICC_MODEL_PATH/Boussinesq/Sphere/Tests/Momentum.hpp )

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

namespace Sphere {

namespace Tests {

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
      #if defined QUICC_SPATIALSCHEME_BLFL || defined QUICC_SPATIALSCHEME_WLFL
         int start = 1;
      #elif defined QUICC_SPATIALSCHEME_BLFM || defined QUICC_SPATIALSCHEME_WLFM
         int start = 0;
      #endif //defined QUICC_SPATIALSCHEME_BLFL || defined QUICC_SPATIALSCHEME_WLFL

      this->defineCoupling(FieldComponents::Spectral::TOR, CouplingInformation::PROGNOSTIC, start, true, false);

      this->defineCoupling(FieldComponents::Spectral::POL, CouplingInformation::PROGNOSTIC, start, true, false);
   }

   void Momentum::setNLComponents()
   {
      this->addNLComponent(FieldComponents::Spectral::TOR, 0);

      this->addNLComponent(FieldComponents::Spectral::POL, 0);
   }

   void Momentum::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id compId) const
   {
      if(compId == FieldComponents::Physical::R)
      {
         rNLComp.setZeros();
      } else if(compId == FieldComponents::Physical::THETA)
      {
         rNLComp.setZeros();
      } else if(compId == FieldComponents::Physical::PHI)
      {
         bool sf = true;

         int nR = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D,Dimensions::Space::PHYSICAL);
         int nTh = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM2D,Dimensions::Space::PHYSICAL);
         int nPh = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM3D,Dimensions::Space::PHYSICAL);

         Array rGrid = Transform::TransformSelector<Dimensions::Transform::TRA1D>::Type::generateGrid(nR);
         Array thGrid = Transform::TransformSelector<Dimensions::Transform::TRA2D>::Type::generateGrid(nTh);
         Array phGrid = Transform::TransformSelector<Dimensions::Transform::TRA3D>::Type::generateGrid(nPh);

         MHDFloat s10t = std::sin(10.0*Math::PI*this->time());
         MHDFloat c10t = std::cos(10.0*Math::PI*this->time());

         MHDFloat funcRct;
         MHDFloat funcRst;
         MHDFloat funcC1Th;
         MHDFloat funcS1Th;
         Array funcPh0 = Array::Ones(nPh);
         //Array funcPh3 = (3.0*phGrid).array().cos() + (3.0*phGrid).array().sin();

         MHDFloat ampct;
         MHDFloat ampst;
         if(sf)
         {
            ampct = (1./1771.)*std::sqrt(1.0/(3.0*Math::PI));
            ampst = (5.*Math::PI/1771.)*std::sqrt(1.0/(3.0*Math::PI));
         } else
         {
            ampct = 2.0*(5./2.)*std::sqrt(3.0/Math::PI);
            ampst = (5./2.)*std::sqrt(3.0/Math::PI);
         }

         MHDFloat r;
         MHDFloat theta;

         rNLComp.rData().setConstant(0);
         nR = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();
         for(int iR = 0; iR < nR; ++iR)
         {
            r = rGrid(this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(iR));
            if(sf)
            {
               funcRct = r*(26565. - 74382.0*std::pow(r,2) - 143451.0*std::pow(r,4) + 57772.0*std::pow(r,6) + 208000.0*std::pow(r,8));
               funcRst = r*(5313.0*std::pow(r,2) - 5313.0*std::pow(r,4) - 5313.0*std::pow(r,6) + 1313.0*std::pow(r,8) + 3200.0*std::pow(r,10));
            } else
            {
               funcRct = r;
               funcRst = r*(-1.0 + std::pow(r,2));
            }

            nTh = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(iR);
            for(int iTh = 0; iTh < nTh; ++iTh)
            {
               theta = thGrid(this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT2D>(iTh, iR));
               //funcC1Th = std::cos(theta);
               funcS1Th = std::sin(theta);

               rNLComp.addProfile(ampct*funcRct*funcS1Th*(100.0+c10t)*funcPh0,iTh,iR);
               rNLComp.addProfile(ampst*funcRst*funcS1Th*s10t*funcPh0,iTh,iR);
            }
         }
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
   }

}
}
}
}
}
