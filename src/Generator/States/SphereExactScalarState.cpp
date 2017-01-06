/** 
 * @file SphereExactScalarState.cpp
 * @brief Source of the implementation of the equation to generate an exact scalar solution in a sphere
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Configuration includes
//

// System includes
//
#include <tr1/cmath>

// External includes
//

// Class include
//
#include "Generator/States/SphereExactScalarState.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "TypeSelectors/TransformSelector.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   SphereExactScalarState::SphereExactScalarState(SharedEquationParameters spEqParams)
      : IScalarEquation(spEqParams), mTypeId(SphereExactStateIds::NOTUSED), mSpecTypeId(SphereExactStateIds::NOTUSED)
   {
   }

   SphereExactScalarState::~SphereExactScalarState()
   {
   }

   void SphereExactScalarState::setIdentity(const PhysicalNames::Id name)
   {
      // Set the name
      this->setName(name);

      // Set the variable requirements
      this->setRequirements();
   }

   void SphereExactScalarState::setStateType(const SphereExactStateIds::Id id)
   {
      this->mTypeId = id;
   }

   void SphereExactScalarState::setSpectralType(const SphereExactStateIds::Id id)
   {
      this->mSpecTypeId = id;
   }

   void SphereExactScalarState::setHarmonicOptions(const SphereExactScalarState::HarmonicModeType& modes)
   {
      this->mSHModes = modes;
   }

   void SphereExactScalarState::setCoupling()
   {
      bool hasNL = false;
      bool hasSource = false;
      if(this->mTypeId != SphereExactStateIds::NOTUSED)
      {
         hasNL = true;
      }

      if(this->mSpecTypeId != SphereExactStateIds::NOTUSED)
      {
         hasSource = true;
      }

      this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::TRIVIAL, 0, hasNL, hasSource, false);
   }

   void SphereExactScalarState::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id compId) const
   {
      // Assert on scalar component is used
      assert(compId == FieldComponents::Physical::SCALAR);

      if(this->mTypeId == SphereExactStateIds::CONSTANT)
      {
         rNLComp.rData().setConstant(42);
      } else if(this->mTypeId == SphereExactStateIds::HARMONIC)
      {
         throw Exception("HARMONIC state is not implemented in physical space");

      } else if(this->mTypeId == SphereExactStateIds::BENCHTEMPC1)
      {
         int nR = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D,Dimensions::Space::PHYSICAL);
         int nTh = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM2D,Dimensions::Space::PHYSICAL);
         int nPh = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM3D,Dimensions::Space::PHYSICAL);

         Array rGrid = Transform::TransformSelector<Dimensions::Transform::TRA1D>::Type::generateGrid(nR);
         Array thGrid = Transform::TransformSelector<Dimensions::Transform::TRA2D>::Type::generateGrid(nTh);
         Array phGrid = Transform::TransformSelector<Dimensions::Transform::TRA3D>::Type::generateGrid(nPh);

         MHDFloat funcR0;
         MHDFloat funcR3;
         MHDFloat funcTh0;
         MHDFloat funcTh3;
         Array funcPh0 = Array::Ones(nPh);
         Array funcPh3 = (3.0*phGrid).array().cos() + (3.0*phGrid).array().sin();

         // Background state is not solved for (no source term but background is imposed explicitly)
         MHDFloat amp0 = 0.0;
         //MHDFloat amp0 = 0.5;
         MHDFloat eps = 1e-5;
         MHDFloat amp3 = (eps/8.0)*std::sqrt(35.0/Math::PI);

         MHDFloat r;
         MHDFloat theta;

         rNLComp.rData().setConstant(0);
         nR = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();
         for(int iR = 0; iR < nR; ++iR)
         {
            r = rGrid(this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(iR));
            funcR0 = (1.0 - std::pow(r,2));
            funcR3 = std::pow(r,3)*(1.0 - std::pow(r,2));

            nTh = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(iR);
            for(int iTh = 0; iTh < nTh; ++iTh)
            {
               theta = thGrid(this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT2D>(iTh, iR));
               funcTh0 = 1.0;
               funcTh3 = std::pow(std::sin(theta),3);

               rNLComp.addProfile(amp0*funcR0*funcTh0*funcPh0,iTh,iR);
               rNLComp.addProfile(amp3*funcR3*funcTh3*funcPh3,iTh,iR);
            }
         }
      } else
      {
         throw Exception("Unknown exact state");
      }
   }

    Datatypes::SpectralScalarType::PointType SphereExactScalarState::sourceTerm(FieldComponents::Spectral::Id compId, const int iN, const int iJ, const int iK) const
    {
      // Assert on scalar component is used
      assert(compId == FieldComponents::Spectral::SCALAR);

      if(this->mSpecTypeId == SphereExactStateIds::HARMONIC)
      {
         #ifdef GEOMHDISCC_SPATIALSCHEME_WLFL
         std::pair<int,int> key = std::make_pair(this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(iK),this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(iJ,iK));
         #endif //GEOMHDISCC_SPATIALSCHEME_WLFL
         #ifdef GEOMHDISCC_SPATIALSCHEME_WLFM
         std::pair<int,int> key = std::make_pair(this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(iJ,iK),this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(iK));
         #endif //GEOMHDISCC_SPATIALSCHEME_WLFM
         if(this->mSHModes.count(key) > 0)
         {
            if(this->mSHModes.find(key)->second.count(iN) > 0)
            {
               return this->mSHModes.find(key)->second.find(iN)->second;
            } else
            {
               return Datatypes::SpectralScalarType::PointType(0.0);
            }
         } else
         {
            return Datatypes::SpectralScalarType::PointType(0);
         }
      } else
      {
         return Datatypes::SpectralScalarType::PointType(0);
      }
    }

   void SphereExactScalarState::setRequirements()
   {
      // Set solver timing
      this->setSolveTiming(SolveTiming::AFTER);

      // Add unknown to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(this->name(), FieldRequirement(true, true, false, false));
   }

}
}
