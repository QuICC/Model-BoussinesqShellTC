/** 
 * @file Transport.cpp
 * @brief Source of the implementation of the transport equation in the Boussinesq rotating thermal convection in a sphere
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
//

// Class include
//
#include MAKE_STR( QUICC_MODEL_PATH/Boussinesq/Sphere/Tests/Transport.hpp )

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "Enums/NonDimensional.hpp"
#include "TypeSelectors/TransformSelector.hpp"

namespace QuICC {

namespace Equations {

namespace Boussinesq {

namespace Sphere {

namespace Tests {

   Transport::Transport(SharedEquationParameters spEqParams)
      : IScalarEquation(spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   Transport::~Transport()
   {
   }

   void Transport::setCoupling()
   {
      this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::PROGNOSTIC, 0, false, false);
   }

   void Transport::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const
   {
      // Assert on scalar component is used
      assert(id == FieldComponents::Physical::SCALAR);

      bool sf = true;

      int nR = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D,Dimensions::Space::PHYSICAL);
      int nTh = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM2D,Dimensions::Space::PHYSICAL);
      int nPh = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM3D,Dimensions::Space::PHYSICAL);

      Array rGrid = Transform::TransformSelector<Dimensions::Transform::TRA1D>::Type::generateGrid(nR);
      Array thGrid = Transform::TransformSelector<Dimensions::Transform::TRA2D>::Type::generateGrid(nTh);
      Array phGrid = Transform::TransformSelector<Dimensions::Transform::TRA3D>::Type::generateGrid(nPh);

      MHDFloat s5t = std::sin(5.0*this->time());
      MHDFloat c5t = std::cos(5.0*this->time());

      MHDFloat funcRct;
      MHDFloat funcRst;
      MHDFloat funcC1Th;
      Array funcPh0 = Array::Ones(nPh);
      //Array funcPh3 = (3.0*phGrid).array().cos() + (3.0*phGrid).array().sin();

      MHDFloat ampct;
      MHDFloat ampst;
      if(sf)
      {
         ampct = (252./84.)*std::sqrt(1.0/(3.0*Math::PI));
         ampst = (5./84.)*std::sqrt(1.0/(3.0*Math::PI));
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
            funcRct = r*(-5.0 + 7.0*std::pow(r,2));
            funcRst = r*(55.0 - 126.0*std::pow(r,2) + 63.0*std::pow(r,4));
         } else
         {
            funcRct = r;
            funcRst = r*(-1.0 + std::pow(r,2));
         }

         nTh = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(iR);
         for(int iTh = 0; iTh < nTh; ++iTh)
         {
            theta = thGrid(this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT2D>(iTh, iR));
            funcC1Th = std::cos(theta);

            rNLComp.addProfile(ampct*funcRct*funcC1Th*(100.0+c5t)*funcPh0,iTh,iR);
            rNLComp.addProfile(ampst*funcRst*funcC1Th*s5t*funcPh0,iTh,iR);
         }
      }
   }

   void Transport::setRequirements()
   {
      // Set temperatur as equation unknown
      this->setName(PhysicalNames::TEMPERATURE);

      // Set solver timing
      this->setSolveTiming(SolveTiming::PROGNOSTIC);

      // Add temperature to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::TEMPERATURE, FieldRequirement(true, true, false, true));

      // Add velocity to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::VELOCITY, FieldRequirement(false, true, true, false, false));
   }

}
}
}
}
}
