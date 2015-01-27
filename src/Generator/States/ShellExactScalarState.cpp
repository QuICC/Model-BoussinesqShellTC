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
#include "TypeSelectors/EquationEigenSelector.hpp"

namespace GeoMHDiSCC {

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
      this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::TRIVIAL, 0, true, false, false, false);
      this->setExplicitTiming(FieldComponents::Spectral::SCALAR, ExplicitTiming::LINEAR);
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
               funcR = std::pow(r,4) + std::pow(r,3) + std::pow(r,2) + r + 1.0;

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
