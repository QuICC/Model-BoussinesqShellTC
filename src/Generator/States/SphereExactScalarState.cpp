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
#include "TypeSelectors/EquationEigenSelector.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   SphereExactScalarState::SphereExactScalarState(SharedEquationParameters spEqParams)
      : IScalarEquation(spEqParams), mTypeId(SphereExactStateIds::CONSTANT)
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

   void SphereExactScalarState::setHarmonicOptions(const std::vector<SphereExactScalarState::HarmonicModeType>& modes)
   {
      this->mSHModes = modes;
   }

   void SphereExactScalarState::setCoupling()
   {
      this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::TRIVIAL, 0, true, false, false);
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
         int nR = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D,Dimensions::Space::PHYSICAL);
         int nTh = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM2D,Dimensions::Space::PHYSICAL);
         int nPh = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM3D,Dimensions::Space::PHYSICAL);

         Array rGrid = Transform::TransformSelector<Dimensions::Transform::TRA1D>::Type::generateGrid(nR);
         Array thGrid = Transform::TransformSelector<Dimensions::Transform::TRA2D>::Type::generateGrid(nTh);
         Array phGrid = Transform::TransformSelector<Dimensions::Transform::TRA3D>::Type::generateGrid(nPh);

         Array sphHarm(nPh);
         MHDFloat funcR;
         typedef std::vector<HarmonicModeType>::const_iterator ModeIt;
         ModeIt it;
         std::pair<ModeIt, ModeIt>  modeRange = std::make_pair(this->mSHModes.begin(), this->mSHModes.end());

         rNLComp.rData().setConstant(0);
         MHDFloat r;
         MHDFloat theta;
         nR = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();
         for(int iR = 0; iR < nR; ++iR)
         {
            r = rGrid(this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(iR));
            nTh = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(iR);
            for(int iTh = 0; iTh < nTh; ++iTh)
            {
               theta = thGrid(this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT2D>(iTh, iR));
               for(it = modeRange.first; it != modeRange.second; ++it)
               {
                  int l = std::tr1::get<0>(*it);
                  int m = std::tr1::get<1>(*it);

                  funcR = 100.0;
                  for(int n = 1; n < 63; ++n)
                  {
                     funcR += std::pow(r,2*n);
                  }
                  funcR *= std::pow(r,l);

                  MHDComplex amplitude = std::tr1::get<2>(*it);

                  // Spherical harmonic Y_l^m
                  sphHarm = SphereExactStateIds::sph_harmonic(amplitude, l, m, theta, phGrid);

                  rNLComp.addProfile(funcR*sphHarm,iTh,iR);
               }
            }
         }
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

         MHDFloat amp0 = 0.5;
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

    Datatypes::SpectralScalarType::PointType SphereExactScalarState::sourceTerm(FieldComponents::Spectral::Id compId, const int i, const int j, const int k) const
    {
      // Assert on scalar component is used
      assert(compId == FieldComponents::Spectral::SCALAR);

      return Datatypes::SpectralScalarType::PointType(0);
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
