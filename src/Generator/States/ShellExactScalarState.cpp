/** 
 * @file ShellExactScalarState.cpp
 * @brief Source of the implementation of the equation to generate an exact scalar solution in a spherical shell
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
#include "Generator/States/ShellExactScalarState.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "TypeSelectors/TransformSelector.hpp"

namespace QuICC {

namespace Equations {

   ShellExactScalarState::ShellExactScalarState(SharedEquationParameters spEqParams)
      : IScalarEquation(spEqParams), mTypeId(ShellExactStateIds::CONSTANT)
   {
   }

   ShellExactScalarState::~ShellExactScalarState()
   {
   }

   void ShellExactScalarState::setIdentity(const PhysicalNames::Id name)
   {
      // Set the name
      this->setName(name);

      // Set the variable requirements
      this->setRequirements();
   }

   void ShellExactScalarState::setStateType(const ShellExactStateIds::Id id)
   {
      this->mTypeId = id;
   }

   void ShellExactScalarState::setHarmonicOptions(const std::vector<ShellExactScalarState::HarmonicModeType>& modes)
   {
      this->mSHModes = modes;
   }

   void ShellExactScalarState::setCoupling()
   {
      this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::TRIVIAL, 0, true, false, false);
   }

   void ShellExactScalarState::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id compId) const
   {
      // Assert on scalar component is used
      assert(compId == FieldComponents::Physical::SCALAR);

      if(this->mTypeId == ShellExactStateIds::CONSTANT)
      {
         rNLComp.rData().setConstant(42);
      } else if(this->mTypeId == ShellExactStateIds::HARMONIC)
      {
         int nR = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D,Dimensions::Space::PHYSICAL);
         int nTh = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM2D,Dimensions::Space::PHYSICAL);
         int nPh = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM3D,Dimensions::Space::PHYSICAL);

         MHDFloat ro = this->eqParams().nd(NonDimensional::RO);
         MHDFloat rratio = this->eqParams().nd(NonDimensional::RRATIO);
         Array rGrid = Transform::TransformSelector<Dimensions::Transform::TRA1D>::Type::generateGrid(nR, ro, rratio);
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
               funcR = 1.0;
               for(int n = 1; n < 10; ++n)
               {
                  funcR += std::pow(r,n);
               }

               theta = thGrid(this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT2D>(iTh, iR));
               for(it = modeRange.first; it != modeRange.second; ++it)
               {
                  int l = std::tr1::get<0>(*it);
                  int m = std::tr1::get<1>(*it);
                  MHDComplex amplitude = std::tr1::get<2>(*it);

                  // Spherical harmonic Y_l^m
                  sphHarm = ShellExactStateIds::sph_harmonic(amplitude, l, m, theta, phGrid);

                  rNLComp.addProfile(funcR*sphHarm,iTh,iR);
               }
            }
         }
      } else if(this->mTypeId == ShellExactStateIds::BENCHTEMPC1)
      {
         int nR = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D,Dimensions::Space::PHYSICAL);
         int nTh = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM2D,Dimensions::Space::PHYSICAL);
         int nPh = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM3D,Dimensions::Space::PHYSICAL);

         MHDFloat ro = this->eqParams().nd(NonDimensional::RO);
         MHDFloat rratio = this->eqParams().nd(NonDimensional::RRATIO);
         Array rGrid = Transform::TransformSelector<Dimensions::Transform::TRA1D>::Type::generateGrid(nR, ro, rratio);
         Array thGrid = Transform::TransformSelector<Dimensions::Transform::TRA2D>::Type::generateGrid(nTh);
         Array phGrid = Transform::TransformSelector<Dimensions::Transform::TRA3D>::Type::generateGrid(nPh);

         MHDFloat funcR;
         MHDFloat funcTh;
         Array funcPh = (4.0*phGrid).array().cos();

         MHDFloat amplitude = 21.0/std::sqrt(17920*Math::PI);

         MHDFloat r;
         MHDFloat x;
         MHDFloat theta;

         rNLComp.rData().setConstant(0);
         nR = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();
         for(int iR = 0; iR < nR; ++iR)
         {
            r = rGrid(this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(iR));
            x = 2.0*r - rratio*ro - ro;
            funcR = 1.0 - 3.0*std::pow(x,2) + 3.0*std::pow(x,4) - std::pow(x,6);

            nTh = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(iR);
            for(int iTh = 0; iTh < nTh; ++iTh)
            {
               theta = thGrid(this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT2D>(iTh, iR));
               funcTh = std::pow(std::sin(theta),4);

               rNLComp.addProfile(amplitude*funcR*funcTh*funcPh,iTh,iR);
            }
         }
      } else
      {
         throw Exception("Unknown exact state");
      }
   }

    Datatypes::SpectralScalarType::PointType ShellExactScalarState::sourceTerm(FieldComponents::Spectral::Id compId, const int i, const int j, const int k) const
    {
      // Assert on scalar component is used
      assert(compId == FieldComponents::Spectral::SCALAR);

      return Datatypes::SpectralScalarType::PointType(0);
    }

   void ShellExactScalarState::setRequirements()
   {
      // Set solver timing
      this->setSolveTiming(SolveTiming::AFTER);

      // Add unknown to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(this->name(), FieldRequirement(true, true, false, false));
   }

}
}
