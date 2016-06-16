/** 
 * @file CartesianExactVectorState.cpp
 * @brief Source of the implementation of the equation to generate an exact scalar solution in cartesian geometries
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "Generator/States/CartesianExactVectorState.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "TypeSelectors/TransformSelector.hpp"
#include "TypeSelectors/EquationEigenSelector.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   CartesianExactVectorState::CartesianExactVectorState(SharedEquationParameters spEqParams)
      : IVectorEquation(spEqParams)
   {
   }

   CartesianExactVectorState::~CartesianExactVectorState()
   {
   }

   void CartesianExactVectorState::setIdentity(const PhysicalNames::Id name)
   {
      // Set the name
      this->setName(name);

      // Set the variable requirements
      this->setRequirements();
   }

   void CartesianExactVectorState::setStateType(const CartesianExactStateIds::Id id)
   {
      if(FieldComponents::Physical::ONE != FieldComponents::Physical::NOTUSED)
      {
         this->mTypeId.insert(std::make_pair(FieldComponents::Physical::ONE, id));
      }

      if(FieldComponents::Physical::TWO != FieldComponents::Physical::NOTUSED)
      {
         this->mTypeId.insert(std::make_pair(FieldComponents::Physical::TWO, id));
      }

      if(FieldComponents::Physical::THREE != FieldComponents::Physical::NOTUSED)
      {
         this->mTypeId.insert(std::make_pair(FieldComponents::Physical::THREE, id));
      }
   }

   void CartesianExactVectorState::setStateType(const FieldComponents::Physical::Id compId, const CartesianExactStateIds::Id id)
   {
      this->mTypeId.insert(std::make_pair(compId, id));
   }

   void CartesianExactVectorState::setModeOptions(const FieldComponents::Physical::Id compId, const MHDFloat a1, const MHDFloat k1, const MHDFloat a2, const MHDFloat k2)
   {
      Array tmp(2);
      tmp(0) = a1;
      tmp(1) = a2;
      this->mModeA.insert(std::make_pair(compId, tmp));

      tmp(0) = k1;
      tmp(1) = k2;
      this->mModeK.insert(std::make_pair(compId, tmp));
   }

   void CartesianExactVectorState::setModeOptions(const FieldComponents::Physical::Id compId, const MHDFloat a1, const MHDFloat k1, const MHDFloat a2, const MHDFloat k2, const MHDFloat a3, const MHDFloat k3)
   {
      Array tmp(3);
      tmp(0) = a1;
      tmp(1) = a2;
      tmp(2) = a3;
      this->mModeA.insert(std::make_pair(compId, tmp));

      tmp(0) = k1;
      tmp(1) = k2;
      tmp(2) = k3;
      this->mModeK.insert(std::make_pair(compId, tmp));
   }

   void CartesianExactVectorState::setCoupling()
   {
      if(FieldComponents::Spectral::ONE != FieldComponents::Spectral::NOTUSED)
      {
         this->defineCoupling(FieldComponents::Spectral::ONE, CouplingInformation::TRIVIAL, 0, true, false, false);
      }

      if(FieldComponents::Spectral::TWO != FieldComponents::Spectral::NOTUSED)
      {
         this->defineCoupling(FieldComponents::Spectral::TWO, CouplingInformation::TRIVIAL, 0, true, false, false);
      }

      if(FieldComponents::Spectral::THREE != FieldComponents::Spectral::NOTUSED)
      {
         this->defineCoupling(FieldComponents::Spectral::THREE, CouplingInformation::TRIVIAL, 0, true, false, false);
      }
   }

   void CartesianExactVectorState::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id compId) const
   {
      CartesianExactStateIds::Id typeId = this->mTypeId.find(compId)->second;

      if(typeId == CartesianExactStateIds::CONSTANT)
      {
         Array modeA = this->mModeA.find(compId)->second;
         rNLComp.rData().setConstant(modeA.prod());

      // Generate divergence free state for toroidal/poloidal decomposition test
      } else if(typeId == CartesianExactStateIds::TORPOLTFF)
      {
         #ifdef GEOMHDISCC_SPATIALDIMENSION_3D
            int nK = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D,Dimensions::Space::PHYSICAL);
            int nJ = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM2D,Dimensions::Space::PHYSICAL);
            int nI = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM3D,Dimensions::Space::PHYSICAL);

            Array gK = Transform::TransformSelector<Dimensions::Transform::TRA1D>::Type::generateGrid(nK);
            Array gJ = Transform::TransformSelector<Dimensions::Transform::TRA2D>::Type::generateGrid(nJ);
            Array gI = Transform::TransformSelector<Dimensions::Transform::TRA3D>::Type::generateGrid(nI);
         #else
            int nK = 1;
            int nJ = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D,Dimensions::Space::PHYSICAL);
            int nI = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM2D,Dimensions::Space::PHYSICAL);

            Array gK = -42*Array::Ones(1);
            Array gJ = Transform::TransformSelector<Dimensions::Transform::TRA1D>::Type::generateGrid(nJ);
            Array gI = Transform::TransformSelector<Dimensions::Transform::TRA2D>::Type::generateGrid(nI);
         #endif //GEOMHDISCC_SPATIALDIMENSION_3D

         MHDFloat k_;
         MHDFloat j_;
         MHDFloat i_;
         nK = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRAND)->dim<Dimensions::Data::DAT3D>();
         for(int iK = 0; iK < nK; ++iK)
         {
            k_ = gK(this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRAND)->idx<Dimensions::Data::DAT3D>(iK));
            nJ = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRAND)->dim<Dimensions::Data::DAT2D>(iK);
            for(int iJ = 0; iJ < nJ; ++iJ)
            {
               j_ = gJ(this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRAND)->idx<Dimensions::Data::DAT2D>(iJ, iK));
               for(int iI = 0; iI < nI; ++iI)
               {
                  i_ = gI(iI);

                  MHDFloat val = 0.0;

                  if(compId == FieldComponents::Physical::X)
                  {
                     //val = -k_*(std::cos(j_) + 4.0*std::sin(j_))*std::sin(i_);
                     val = -4.0*k_*std::sin(4.0*i_);
                     rNLComp.setPoint(val, iI, iJ, iK);
                  } else if(compId == FieldComponents::Physical::Y)
                  {
                     //val = k_*(4.0*std::cos(j_) + std::sin(j_))*std::cos(i_);
                     val = 4.0*k_*std::cos(i_);
                     rNLComp.setPoint(val, iI, iJ, iK);
                  } else if(compId == FieldComponents::Physical::Z)
                  {
                     for(int sI = 0; sI < 5; sI++)
                     {
                        for(int sJ = 0; sJ < 5; sJ++)
                        {
                           val += -(sI*sI + sJ*sJ)*(std::cos(sI*i_)+std::sin(sI*i_))*(std::cos(sJ*j_)+std::sin(sJ*j_));
                        }
                     }
                     rNLComp.setPoint(val, iI, iJ, iK);
                  }
               }
            }
         }

      } else
      {
         Array modeA = this->mModeA.find(compId)->second;
         Array modeK = this->mModeK.find(compId)->second;
         #ifdef GEOMHDISCC_SPATIALDIMENSION_3D
            int nK = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D,Dimensions::Space::PHYSICAL);
            int nJ = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM2D,Dimensions::Space::PHYSICAL);
            int nI = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM3D,Dimensions::Space::PHYSICAL);

            Array gK = Transform::TransformSelector<Dimensions::Transform::TRA1D>::Type::generateGrid(nK);
            Array gJ = Transform::TransformSelector<Dimensions::Transform::TRA2D>::Type::generateGrid(nJ);
            Array gI = Transform::TransformSelector<Dimensions::Transform::TRA3D>::Type::generateGrid(nI);
         #else
            int nK = 1;
            int nJ = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D,Dimensions::Space::PHYSICAL);
            int nI = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM2D,Dimensions::Space::PHYSICAL);

            Array gK = -42*Array::Ones(1);
            Array gJ = Transform::TransformSelector<Dimensions::Transform::TRA1D>::Type::generateGrid(nJ);
            Array gI = Transform::TransformSelector<Dimensions::Transform::TRA2D>::Type::generateGrid(nI);
         #endif //GEOMHDISCC_SPATIALDIMENSION_3D

         MHDFloat k_;
         MHDFloat j_;
         MHDFloat i_;
         nK = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRAND)->dim<Dimensions::Data::DAT3D>();
         for(int iK = 0; iK < nK; ++iK)
         {
            k_ = gK(this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRAND)->idx<Dimensions::Data::DAT3D>(iK));
            nJ = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRAND)->dim<Dimensions::Data::DAT2D>(iK);
            for(int iJ = 0; iJ < nJ; ++iJ)
            {
               j_ = gJ(this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRAND)->idx<Dimensions::Data::DAT2D>(iJ, iK));
               for(int iI = 0; iI < nI; ++iI)
               {
                  i_ = gI(iI);

                  MHDFloat val = 0.0;
                  if(static_cast<int>(typeId) < 100)
                  {
                     Array grid(3);
                     grid(0) = k_;
                     grid(1) = j_;
                     grid(2) = i_;
                     val = CartesianExactStateIds::exact3D(typeId, modeA, modeK, grid);
                  } else if(static_cast<int>(typeId) >= 100)
                  {
                     Array grid(2);
                     grid(0) = j_;
                     grid(1) = i_;
                     val = CartesianExactStateIds::exact2D(typeId, modeA, modeK, grid);
                  }

                  rNLComp.setPoint(val, iI, iJ, iK);
               }
            }
         }
      }
   }

    Datatypes::SpectralScalarType::PointType CartesianExactVectorState::sourceTerm(FieldComponents::Spectral::Id compId, const int i, const int j, const int k) const
    {
      return Datatypes::SpectralScalarType::PointType(0);
    }

   void CartesianExactVectorState::setRequirements()
   {
      // Set solver timing
      this->setSolveTiming(SolveTiming::AFTER);

      // Add unknown to requirements: is scalar?, need spectral?, need physical?, need diff?(, need curl?)
      this->mRequirements.addField(this->name(), FieldRequirement(false, true, true, false));
   }

   void CartesianExactVectorState::setNLComponents()
   {
      if(FieldComponents::Spectral::ONE != FieldComponents::Spectral::NOTUSED)
      {
         this->addNLComponent(FieldComponents::Spectral::ONE, 0);
      }

      if(FieldComponents::Spectral::TWO != FieldComponents::Spectral::NOTUSED)
      {
         this->addNLComponent(FieldComponents::Spectral::TWO, 0);
      }

      if(FieldComponents::Spectral::THREE != FieldComponents::Spectral::NOTUSED)
      {
         this->addNLComponent(FieldComponents::Spectral::THREE, 0);
      }
   }

}
}
