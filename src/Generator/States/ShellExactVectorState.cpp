/** 
 * @file ShellExactVectorState.cpp
 * @brief Source of the implementation of the equation to generate an exact vector solution in a shell
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
#include "Generator/States/ShellExactVectorState.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "TypeSelectors/TransformSelector.hpp"
#include "TypeSelectors/EquationEigenSelector.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   ShellExactVectorState::ShellExactVectorState(SharedEquationParameters spEqParams)
      : IVectorEquation(spEqParams), mTypeId(CONSTANT)
   {
   }

   ShellExactVectorState::~ShellExactVectorState()
   {
   }

   void ShellExactVectorState::setIdentity(const PhysicalNames::Id name)
   {
      // Set the name
      this->setName(name);

      // Set toroidal and poloidal components
      this->mSpectralIds.push_back(FieldComponents::Spectral::ONE);
      this->mSpectralIds.push_back(FieldComponents::Spectral::TWO);

      // Set the variable requirements
      this->setRequirements();
   }

   void ShellExactVectorState::setStateType(const ShellExactVectorState::StateTypeId id)
   {
      this->mTypeId = id;
   }

   void ShellExactVectorState::setHarmonicOptions(const std::vector<std::tr1::tuple<int,int,MHDComplex> >& rModes, const std::vector<std::tr1::tuple<int,int,MHDComplex> >& tModes, const std::vector<std::tr1::tuple<int,int,MHDComplex> >& pModes)
   {
      this->mRSHModes = rModes;
      this->mTSHModes = tModes;
      this->mPSHModes = pModes;
   }

   void ShellExactVectorState::setCoupling()
   {
      SpectralComponent_range specRange = this->spectralRange();
      SpectralComponent_iterator specIt;
      for(specIt = specRange.first; specIt != specRange.second; ++specIt)
      {
         this->defineCoupling(*specIt, CouplingInformation::TRIVIAL, 0, true, false, false);
      }
   }

   void ShellExactVectorState::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id compId) const
   {
      if(this->mTypeId == CONSTANT)
      {
         rNLComp.rData().setConstant(42);
      } else if(this->mTypeId == HARMONIC)
      {
         int nR = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D,Dimensions::Space::PHYSICAL);
         int nTh = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM2D,Dimensions::Space::PHYSICAL);
         int nPh = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM3D,Dimensions::Space::PHYSICAL);

         Array rGrid = Transform::TransformSelector<Dimensions::Transform::TRA1D>::Type::generateGrid(nR, 1, 0.35);
         Array thGrid = Transform::TransformSelector<Dimensions::Transform::TRA2D>::Type::generateGrid(nTh);
         Array phGrid = Transform::TransformSelector<Dimensions::Transform::TRA3D>::Type::generateGrid(nPh);

         Array funcPh(nPh);
         MHDFloat funcR(nR);
         MHDFloat funcTh;
         typedef std::vector<std::tr1::tuple<int,int,MHDComplex> >::const_iterator ModeIt;
         std::pair<ModeIt, ModeIt>  modeRange;
         if(compId == FieldComponents::Physical::ONE)
         {
            modeRange.first = this->mRSHModes.begin();
            modeRange.second = this->mRSHModes.end();
         } else if(compId == FieldComponents::Physical::TWO)
         {
            modeRange.first = this->mTSHModes.begin();
            modeRange.second = this->mTSHModes.end();
         } else
         {
            modeRange.first = this->mPSHModes.begin();
            modeRange.second = this->mPSHModes.end();
         }
         rNLComp.rData().setConstant(0);
         for(int iR = 0; iR < nR; ++iR)
         {
            for(int iTh = 0; iTh < nTh; ++iTh)
            {
               funcR = 1.0;

               for(ModeIt it = modeRange.first; it != modeRange.second; ++it)
               {
                  int l = std::tr1::get<0>(*it);
                  int m = std::tr1::get<1>(*it);
                  MHDFloat re = std::tr1::get<2>(*it).real();
                  MHDFloat im = std::tr1::get<2>(*it).imag();

                  // Spherical harmonic Y_l^m
                  funcPh = re*(static_cast<MHDFloat>(m)*phGrid).array().cos() + im*(static_cast<MHDFloat>(m)*phGrid).array().sin();
                  funcTh = std::tr1::sph_legendre(l,m, thGrid(iTh));
                  rNLComp.addProfile(funcPh*funcR*funcTh,iTh,iR);
               }
            }
         }
      } else
      {
         throw Exception("Unknown exact state");
      }
   }

    Datatypes::SpectralScalarType::PointType ShellExactVectorState::sourceTerm(FieldComponents::Spectral::Id compId, const int i, const int j, const int k) const
    {
      return Datatypes::SpectralScalarType::PointType(0);
    }

   void ShellExactVectorState::setRequirements()
   {
      // Add unknown to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(this->name(), FieldRequirement(false, true, false, false));
   }

}
}
