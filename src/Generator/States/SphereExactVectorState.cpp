/** 
 * @file SphereExactVectorState.cpp
 * @brief Source of the implementation of the equation to generate an exact vector solution in a sphere
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
#include "Generator/States/SphereExactVectorState.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "TypeSelectors/TransformSelector.hpp"

namespace QuICC {

namespace Equations {

   SphereExactVectorState::SphereExactVectorState(SharedEquationParameters spEqParams)
      : IVectorEquation(spEqParams), mTypeId(SphereExactStateIds::NOTUSED), mSpecTypeId(SphereExactStateIds::NOTUSED)
   {
   }

   SphereExactVectorState::~SphereExactVectorState()
   {
   }

   void SphereExactVectorState::setIdentity(const PhysicalNames::Id name)
   {
      // Set the name
      this->setName(name);

      // Set the variable requirements
      this->setRequirements();
   }

   void SphereExactVectorState::setStateType(const SphereExactStateIds::Id id)
   {
      this->mTypeId = id;
   }

   void SphereExactVectorState::setSpectralType(const SphereExactStateIds::Id id)
   {
      this->mSpecTypeId = id;
   }

   void SphereExactVectorState::setHarmonicOptions(const FieldComponents::Spectral::Id compId, const SphereExactVectorState::HarmonicModeType& modes)
   {
      this->mSHModes.insert(std::make_pair(compId, modes));
   }

   void SphereExactVectorState::setCoupling()
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

      if(FieldComponents::Spectral::ONE != FieldComponents::Spectral::NOTUSED)
      {
         this->defineCoupling(FieldComponents::Spectral::ONE, CouplingInformation::TRIVIAL, 0, hasNL, hasSource, false, false);
      }

      if(FieldComponents::Spectral::TWO != FieldComponents::Spectral::NOTUSED)
      {
         this->defineCoupling(FieldComponents::Spectral::TWO, CouplingInformation::TRIVIAL, 0, hasNL, hasSource, false, false);
      }

      if(FieldComponents::Spectral::THREE != FieldComponents::Spectral::NOTUSED)
      {
         this->defineCoupling(FieldComponents::Spectral::THREE, CouplingInformation::TRIVIAL, 0, hasNL, hasSource, false, false);
      }
   }

   void SphereExactVectorState::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id compId) const
   {
      SphereExactStateIds::Id typeId = this->mTypeId;

      // Initialize to zero
      rNLComp.rData().setConstant(0);

      if(typeId == SphereExactStateIds::CONSTANT)
      {
         rNLComp.rData().setConstant(42);
      } else if(typeId == SphereExactStateIds::HARMONIC)
      {
         throw Exception("HARMONIC state is not implemented in physical space");

      } else if(typeId == SphereExactStateIds::BENCHVELC1)
      {
         rNLComp.rData().setConstant(1e-16);

      } else if(typeId == SphereExactStateIds::BENCHVELC2)
      {
         int nR = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D,Dimensions::Space::PHYSICAL);
         int nTh = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM2D,Dimensions::Space::PHYSICAL);
         int nPh = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM3D,Dimensions::Space::PHYSICAL);

         Array rGrid = Transform::TransformSelector<Dimensions::Transform::TRA1D>::Type::generateGrid(nR);
         Array thGrid = Transform::TransformSelector<Dimensions::Transform::TRA2D>::Type::generateGrid(nTh);
         Array phGrid = Transform::TransformSelector<Dimensions::Transform::TRA3D>::Type::generateGrid(nPh);

         Array funcPh = Array::Ones(nPh);
         Array funcPhC1 = (phGrid).array().cos();
         Array funcPhS1 = (phGrid).array().sin();
         MHDFloat funcR = 1.0;
         MHDFloat amplitude = 1.0;
         MHDFloat funcTh = 1.0;

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

               if(compId == FieldComponents::Physical::R)
               {
                  amplitude = 0.0;
                  rNLComp.addProfile(amplitude*funcPh,iTh,iR);
               } else if(compId == FieldComponents::Physical::THETA)
               {
                  amplitude = -10.0/(7.0*std::sqrt(3.0));
                  funcTh = std::cos(theta);

                  funcR = 3.0*std::pow(r,2)*(-147.0 + 343.0*std::pow(r,2) - 217.0*std::pow(r,4) + 29.0*std::pow(r,6));
                  rNLComp.addProfile(amplitude*funcR*funcTh*funcPhC1,iTh,iR);

                  funcR = 14.0*std::pow(r,2)*(-9.0 - 125.0*std::pow(r,2) + 39.0*std::pow(r,4) + 27.0*std::pow(r,6));
                  rNLComp.addProfile(amplitude*funcR*funcTh*funcPhS1,iTh,iR);
               } else if(compId == FieldComponents::Physical::PHI)
               {
                  amplitude = -5.0/5544.;

                  funcR = 7.0*r*(43700.0 - 58113.0*std::pow(r,2) - 15345.0*std::pow(r,4) + 1881.0*std::pow(r,6) + 20790.0*std::pow(r,8));
                  funcTh = std::sin(theta);
                  rNLComp.addProfile(amplitude*funcR*funcTh*funcPh,iTh,iR);

                  funcR = 7.0*1485*std::pow(r,3)*(-9.0 + 115.0*std::pow(r,2) - 167.0*std::pow(r,4) + 70.0*std::pow(r,6));
                  funcTh = std::sin(3.0*theta);
                  rNLComp.addProfile(amplitude*funcR*funcTh*funcPh,iTh,iR);

                  funcR = 528*std::sqrt(3)*std::pow(r,2)*14*(-9.0 - 125.0*std::pow(r,2) + 39.0*std::pow(r,4) + 27.0*std::pow(r,6));
                  funcTh = std::cos(2.0*theta);
                  rNLComp.addProfile(amplitude*funcR*funcTh*funcPhC1,iTh,iR);

                  funcR = 528*std::sqrt(3)*std::pow(r,2)*3*(147.0 - 343.0*std::pow(r,2) + 217.0*std::pow(r,4) - 29.0*std::pow(r,6));
                  funcTh = std::cos(2.0*theta);
                  rNLComp.addProfile(amplitude*funcR*funcTh*funcPhS1,iTh,iR);
               }
            }
         }

      } else if(typeId == SphereExactStateIds::BENCHMAGC2)
      {
         int nR = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D,Dimensions::Space::PHYSICAL);
         int nTh = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM2D,Dimensions::Space::PHYSICAL);
         int nPh = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM3D,Dimensions::Space::PHYSICAL);

         Array rGrid = Transform::TransformSelector<Dimensions::Transform::TRA1D>::Type::generateGrid(nR);
         Array thGrid = Transform::TransformSelector<Dimensions::Transform::TRA2D>::Type::generateGrid(nTh);
         Array phGrid = Transform::TransformSelector<Dimensions::Transform::TRA3D>::Type::generateGrid(nPh);

         Array funcPh = Array::Ones(nPh);
         Array funcPhCS = (phGrid).array().cos() + (phGrid).array().sin();
         Array funcPhC_S = (phGrid).array().cos() - (phGrid).array().sin();
         MHDFloat funcR = 1.0;
         MHDFloat funcTh = 1.0;
         MHDFloat amplitude = 1.0;

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

               if(compId == FieldComponents::Physical::R)
               {
                  amplitude = 0.0;
                  rNLComp.addProfile(amplitude*funcR*funcTh*funcPh,iTh,iR);
               } else if(compId == FieldComponents::Physical::THETA)
               {
                  amplitude = -3.0/2.0;
                  funcR = r*(-1.0 + 4.0*std::pow(r,2) - 6.0*std::pow(r,4) + 3.0*std::pow(r,6));
                  funcTh = 1.0;
                  rNLComp.addProfile(amplitude*funcR*funcTh*funcPhCS,iTh,iR);
               } else if(compId == FieldComponents::Physical::PHI)
               {
                  amplitude = -3.0/4.0;

                  funcR = 3.0*std::pow(r,2)*(-1.0 + std::pow(r,2))*(2.0 - 5.0*std::pow(r,2) + 4.0*std::pow(r,4));
                  funcTh = std::cos(theta)*std::sin(theta);
                  rNLComp.addProfile(amplitude*funcR*funcTh*funcPh,iTh,iR);

                  funcR = 2.0*r*(-1.0 + std::pow(r,2))*(1.0 - 3.0*std::pow(r,2) + 3.0*std::pow(r,4));
                  funcTh = std::cos(theta);
                  rNLComp.addProfile(amplitude*funcR*funcTh*funcPhC_S,iTh,iR);
               }

            }
         }
      } else
      {
         throw Exception("Unknown exact state");
      }
   }

    Datatypes::SpectralScalarType::PointType SphereExactVectorState::sourceTerm(FieldComponents::Spectral::Id compId, const int iN, const int iJ, const int iK) const
    {
      if(this->mSpecTypeId == SphereExactStateIds::HARMONIC && this->mSHModes.count(compId) > 0)
      {
         #ifdef QUICC_SPATIALSCHEME_WLFL
         std::pair<int,int> key = std::make_pair(this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(iK),this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(iJ,iK));
         #endif //QUICC_SPATIALSCHEME_WLFL
         #ifdef QUICC_SPATIALSCHEME_WLFM
         std::pair<int,int> key = std::make_pair(this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(iJ,iK),this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(iK));
         #endif //QUICC_SPATIALSCHEME_WLFM
         if(this->mSHModes.find(compId)->second.count(key) > 0)
         {
            if(this->mSHModes.find(compId)->second.find(key)->second.count(iN) > 0)
            {
               return this->mSHModes.find(compId)->second.find(key)->second.find(iN)->second;
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

   void SphereExactVectorState::setRequirements()
   {
      // Set solver timing
      this->setSolveTiming(SolveTiming::AFTER);

      // Add unknown to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(this->name(), FieldRequirement(false, true, false, false));
   }

   void SphereExactVectorState::setNLComponents()
   {
      if(FieldComponents::Spectral::ONE != FieldComponents::Spectral::NOTUSED)
      {
         this->addNLComponent(FieldComponents::Spectral::ONE, 0);
      }

      if(FieldComponents::Spectral::TWO != FieldComponents::Spectral::NOTUSED)
      {
         this->addNLCom
         onent(FieldComponents::Spectral::TWO, 0);
      }

      if(FieldComponents::Spectral::THREE != FieldComponents::Spectral::NOTUSED)
      {
         this->addNLComponent(FieldComponents::Spectral::THREE, 0);
      }
   }

}
}
