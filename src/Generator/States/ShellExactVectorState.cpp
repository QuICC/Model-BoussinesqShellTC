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
      : IVectorEquation(spEqParams)
   {
   }

   ShellExactVectorState::~ShellExactVectorState()
   {
   }

   void ShellExactVectorState::setIdentity(const PhysicalNames::Id name)
   {
      // Set the name
      this->setName(name);

      // Set the variable requirements
      this->setRequirements();
   }

   void ShellExactVectorState::setStateType(const FieldComponents::Physical::Id compId, const ShellExactStateIds::Id id)
   {
      this->mTypeId.insert(std::make_pair(compId, id));
   }

   void ShellExactVectorState::setHarmonicOptions(const FieldComponents::Physical::Id compId, const std::vector<ShellExactVectorState::HarmonicModeType>& modes)
   {
      this->mSHModes.insert(std::make_pair(compId, modes));
   }

   void ShellExactVectorState::setCoupling()
   {
      if(FieldComponents::Spectral::ONE != FieldComponents::Spectral::NOTUSED)
      {
         this->defineCoupling(FieldComponents::Spectral::ONE, CouplingInformation::TRIVIAL, 0, true, false, false, false);
      }

      if(FieldComponents::Spectral::TWO != FieldComponents::Spectral::NOTUSED)
      {
         this->defineCoupling(FieldComponents::Spectral::TWO, CouplingInformation::TRIVIAL, 0, true, false, false, false);
      }

      if(FieldComponents::Spectral::THREE != FieldComponents::Spectral::NOTUSED)
      {
         this->defineCoupling(FieldComponents::Spectral::THREE, CouplingInformation::TRIVIAL, 0, true, false, false, false);
      }
   }

   void ShellExactVectorState::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id compId) const
   {
      ShellExactStateIds::Id typeId = this->mTypeId.find(compId)->second;

      if(typeId == ShellExactStateIds::CONSTANT)
      {
         rNLComp.rData().setConstant(42);
      } else if(typeId == ShellExactStateIds::HARMONIC)
      {
         int nR = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D,Dimensions::Space::PHYSICAL);
         int nTh = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM2D,Dimensions::Space::PHYSICAL);
         int nPh = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM3D,Dimensions::Space::PHYSICAL);

         Array rGrid = Transform::TransformSelector<Dimensions::Transform::TRA1D>::Type::generateGrid(nR, 1, 0.35);
         Array thGrid = Transform::TransformSelector<Dimensions::Transform::TRA2D>::Type::generateGrid(nTh);
         Array phGrid = Transform::TransformSelector<Dimensions::Transform::TRA3D>::Type::generateGrid(nPh);

         Array sphHarm(nPh);
         MHDFloat funcR;
         typedef std::vector<HarmonicModeType>::const_iterator ModeIt;
         std::pair<ModeIt, ModeIt>  modeRange;

         modeRange.first = this->mSHModes.find(compId)->second.begin();
         modeRange.second = this->mSHModes.find(compId)->second.end();

         rNLComp.rData().setConstant(0);
         for(int iR = 0; iR < nR; ++iR)
         {
            funcR = 1.0;
            funcR = rGrid(iR);

            for(int iTh = 0; iTh < nTh; ++iTh)
            {
               for(ModeIt it = modeRange.first; it != modeRange.second; ++it)
               {
                  int l = std::tr1::get<0>(*it);
                  int m = std::tr1::get<1>(*it);
                  MHDComplex amplitude = std::tr1::get<2>(*it);

                  // Spherical harmonic Y_l^m
                  sphHarm = ShellExactStateIds::sph_harmonic(amplitude, l, m, thGrid(iTh), phGrid);
                  
                  rNLComp.addProfile(funcR*sphHarm,iTh,iR);
               }
            }
         }

      } else if(typeId == ShellExactStateIds::TORPOLT11P11)
      {
         int nR = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D,Dimensions::Space::PHYSICAL);
         int nTh = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM2D,Dimensions::Space::PHYSICAL);
         int nPh = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM3D,Dimensions::Space::PHYSICAL);

         Array rGrid = Transform::TransformSelector<Dimensions::Transform::TRA1D>::Type::generateGrid(nR, 1, 0.35);
         Array thGrid = Transform::TransformSelector<Dimensions::Transform::TRA2D>::Type::generateGrid(nTh);
         Array phGrid = Transform::TransformSelector<Dimensions::Transform::TRA3D>::Type::generateGrid(nPh);

         Array funcPh(nPh);
         MHDFloat funcR = 1.0;
         MHDFloat funcTh = 1.0;

         rNLComp.rData().setConstant(0);
         for(int iR = 0; iR < nR; ++iR)
         {
            for(int iTh = 0; iTh < nTh; ++iTh)
            {
               // Toroidal part
               if(compId == FieldComponents::Physical::R)
               {
                  funcR = 0.0;
                  funcTh = 0.0;
                  funcPh.setConstant(0.0);

               } else if(compId == FieldComponents::Physical::THETA)
               {
                  funcR = std::sqrt(3.0/(2.0*Math::PI));
                  funcTh = 1.0;
                  funcPh = phGrid.array().cos() + phGrid.array().sin();

               } else if(compId == FieldComponents::Physical::PHI)
               {
                  funcR = std::sqrt(3.0/(2.0*Math::PI));
                  funcTh = std::cos(thGrid(iTh));
                  funcPh = phGrid.array().cos() - phGrid.array().sin();
               }
               rNLComp.addProfile(funcR*funcTh*funcPh, iTh, iR);

               // Poloidal part
               if(compId == FieldComponents::Physical::R)
               {
                  funcR = std::sqrt(6.0/Math::PI)/rGrid(iR);
                  funcTh = std::sin(thGrid(iTh));
                  funcPh = -phGrid.array().cos() + phGrid.array().sin();

               } else if(compId == FieldComponents::Physical::THETA)
               {
                  funcR = std::sqrt(3.0/(2.0*Math::PI))/rGrid(iR);
                  funcTh = std::cos(thGrid(iTh));
                  funcPh = -phGrid.array().cos() + phGrid.array().sin();

               } else if(compId == FieldComponents::Physical::PHI)
               {
                  funcR = std::sqrt(3.0/(2.0*Math::PI))/rGrid(iR);
                  funcTh = 1.0;
                  funcPh = phGrid.array().cos() + phGrid.array().sin();
               }
               rNLComp.addProfile(funcR*funcTh*funcPh, iTh, iR);
            }
         }

      } else if(typeId == ShellExactStateIds::TORPOLT54P43)
      {
         int nR = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D,Dimensions::Space::PHYSICAL);
         int nTh = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM2D,Dimensions::Space::PHYSICAL);
         int nPh = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM3D,Dimensions::Space::PHYSICAL);

         Array rGrid = Transform::TransformSelector<Dimensions::Transform::TRA1D>::Type::generateGrid(nR, 1, 0.35);
         Array thGrid = Transform::TransformSelector<Dimensions::Transform::TRA2D>::Type::generateGrid(nTh);
         Array phGrid = Transform::TransformSelector<Dimensions::Transform::TRA3D>::Type::generateGrid(nPh);

         Array funcPh(nPh);
         MHDFloat funcR = 1.0;
         MHDFloat funcTh = 1.0;

         rNLComp.rData().setConstant(0);
         for(int iR = 0; iR < nR; ++iR)
         {
            for(int iTh = 0; iTh < nTh; ++iTh)
            {
               // Toroidal part
               if(compId == FieldComponents::Physical::R)
               {
                  funcR = 0.0;
                  funcTh = 0.0;
                  funcPh.setConstant(0.0);

               } else if(compId == FieldComponents::Physical::THETA)
               {
                  funcR = -(3.0/2.0)*std::sqrt(385.0/(2.0*Math::PI));
                  funcTh = std::cos(thGrid(iTh))*std::pow(std::sin(thGrid(iTh)),3);
                  funcPh = (4.0*phGrid).array().cos() + (4.0*phGrid).array().sin();

               } else if(compId == FieldComponents::Physical::PHI)
               {
                  funcR = -(3.0/16.0)*std::sqrt(385.0/(2.0*Math::PI));
                  funcTh = (3.0 + 5.0*std::cos(2.0*thGrid(iTh)))*std::pow(std::sin(thGrid(iTh)),3);
                  funcPh = (4.0*phGrid).array().cos() - (4.0*phGrid).array().sin();
               }
               rNLComp.addProfile(funcR*funcTh*funcPh, iTh, iR);

               // Poloidal part
               if(compId == FieldComponents::Physical::R)
               {
                  funcR = -15.0*std::sqrt(35.0/Math::PI)*std::pow(rGrid(iR),2)*(1.0 + rGrid(iR));
                  funcTh = std::cos(thGrid(iTh))*std::pow(std::sin(thGrid(iTh)),3);
                  funcPh = (3.0*phGrid).array().cos() - (3.0*phGrid).array().sin();

               } else if(compId == FieldComponents::Physical::THETA)
               {
                  funcR = -(3.0/4.0)*std::sqrt(35.0/Math::PI)*std::pow(rGrid(iR),2)*(4.0 + 5.0*rGrid(iR));
                  funcTh = (1.0 + 2.0*std::cos(2.0*thGrid(iTh)))*std::pow(std::sin(thGrid(iTh)),2);
                  funcPh = (3.0*phGrid).array().cos() - (3.0*phGrid).array().sin();

               } else if(compId == FieldComponents::Physical::PHI)
               {
                  funcR = (9.0/4.0)*std::sqrt(35.0/Math::PI)*std::pow(rGrid(iR),2)*(4.0 + 5.0*rGrid(iR));
                  funcTh = std::cos(thGrid(iTh))*std::pow(std::sin(thGrid(iTh)),2);
                  funcPh = (3.0*phGrid).array().cos() + (3.0*phGrid).array().sin();
               }
               rNLComp.addProfile(funcR*funcTh*funcPh, iTh, iR);
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
      // Set solver timing
      this->setSolveTiming(SolveTiming::AFTER);

      // Add unknown to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(this->name(), FieldRequirement(false, true, false, false));
   }

}
}
