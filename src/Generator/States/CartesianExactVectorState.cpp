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

   void CartesianExactVectorState::setStateType(const FieldComponents::Physical::Id compId, const CartesianExactStateIds::Id id)
   {
      this->mTypeId.insert(std::make_pair(compId, id));
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

   void CartesianExactVectorState::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id compId) const
   {
      CartesianExactStateIds::Id typeId = this->mTypeId.find(compId)->second;
      Array modeA = this->mModeA.find(compId)->second;
      Array modeK = this->mModeK.find(compId)->second;

      if(typeId == CartesianExactStateIds::CONSTANT)
      {
         rNLComp.rData().setConstant(modeA(0)*modeA(1)*modeA(2));
      } else
      {
         int nK = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D,Dimensions::Space::PHYSICAL);
         int nJ = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM2D,Dimensions::Space::PHYSICAL);
         int nI = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM3D,Dimensions::Space::PHYSICAL);

         Array gK = Transform::TransformSelector<Dimensions::Transform::TRA1D>::Type::generateGrid(nK);
         Array gJ = Transform::TransformSelector<Dimensions::Transform::TRA2D>::Type::generateGrid(nJ);
         Array gI = Transform::TransformSelector<Dimensions::Transform::TRA3D>::Type::generateGrid(nI);

         MHDFloat k_;
         MHDFloat j_;
         MHDFloat i_;
         nK = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();
         for(int iK = 0; iK < nK; ++iK)
         {
            k_ = gK(this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(iK));
            nJ = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(iK);
            for(int iJ = 0; iJ < nJ; ++iJ)
            {
               j_ = gJ(this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT2D>(iJ, iK));
               for(int iI = 0; iI < nI; ++iI)
               {
                  i_ = gI(iI);

                  MHDFloat valJ = 0.0;
                  MHDFloat valI = 0.0;
                  MHDFloat valK = 0.0;

                  // Generate solutions for TTT geometries
                  if(static_cast<int>(typeId) > 9 && static_cast<int>(typeId) < 20)
                  {
                     if(typeId == CartesianExactStateIds::POLYPOLYPOLY)
                     {
                        valK = CartesianExactStateIds::poly(modeA(0),modeK(0),k_);
                        valJ = CartesianExactStateIds::poly(modeA(1),modeK(1),j_);
                        valI = CartesianExactStateIds::poly(modeA(2),modeK(2),i_);

                     } else
                     {
                        throw Exception("Unknown exact state");
                     }

                     // Generate solutions for TFT geometries
                  } else if(static_cast<int>(typeId) > 19 && static_cast<int>(typeId) < 30)
                  {
                     if(typeId == CartesianExactStateIds::POLYCOSPOLY)
                     {
                        valK = CartesianExactStateIds::poly(modeA(0),modeK(0),k_);
                        valJ = CartesianExactStateIds::cos(modeA(1),modeK(1),j_);
                        valI = CartesianExactStateIds::poly(modeA(2),modeK(2),i_);

                     } else if(typeId == CartesianExactStateIds::POLYSINPOLY)
                     {
                        valK = CartesianExactStateIds::poly(modeA(0),modeK(0),k_);
                        valJ = CartesianExactStateIds::sin(modeA(1),modeK(1),j_);
                        valI = CartesianExactStateIds::poly(modeA(2),modeK(2),i_);

                     } else
                     {
                        throw Exception("Unknown exact state");
                     }

                     // Generate solutions for TFF geometries
                  } else if(static_cast<int>(typeId) > 29 && static_cast<int>(typeId) < 50)
                  {
                     if(typeId == CartesianExactStateIds::POLYCOSCOS)
                     {
                        valK = CartesianExactStateIds::poly(modeA(0),modeK(0),k_);
                        valJ = CartesianExactStateIds::cos(modeA(1),modeK(1),j_);
                        valI = CartesianExactStateIds::cos(modeA(2),modeK(2),i_);

                     } else if(typeId == CartesianExactStateIds::POLYSINSIN)
                     {
                        valK = CartesianExactStateIds::poly(modeA(0),modeK(0),k_);
                        valJ = CartesianExactStateIds::sin(modeA(1),modeK(1),j_);
                        valI = CartesianExactStateIds::sin(modeA(2),modeK(2),i_);

                     } else if(typeId == CartesianExactStateIds::POLYSINCOS)
                     {
                        valK = CartesianExactStateIds::poly(modeA(0),modeK(0),k_);
                        valJ = CartesianExactStateIds::sin(modeA(1),modeK(1),j_);
                        valI = CartesianExactStateIds::cos(modeA(2),modeK(2),i_);

                     } else if(typeId == CartesianExactStateIds::POLYCOSSIN)
                     {
                        valK = CartesianExactStateIds::poly(modeA(0),modeK(0),k_);
                        valJ = CartesianExactStateIds::cos(modeA(1),modeK(1),j_);
                        valI = CartesianExactStateIds::sin(modeA(2),modeK(2),i_);

                     } else
                     {
                        throw Exception("Unknown exact state");
                     }

                     // Generate solutions for FFF geometries
                  } else if(static_cast<int>(typeId) > 39 && static_cast<int>(typeId) < 50)
                  {
                     if(typeId == CartesianExactStateIds::COSCOSCOS)
                     {
                        valK = CartesianExactStateIds::cos(modeA(0),modeK(0),k_);
                        valJ = CartesianExactStateIds::cos(modeA(1),modeK(1),j_);
                        valI = CartesianExactStateIds::cos(modeA(2),modeK(2),i_);

                     } else if(typeId == CartesianExactStateIds::SINSINSIN)
                     {
                        valK = CartesianExactStateIds::sin(modeA(0),modeK(0),k_);
                        valJ = CartesianExactStateIds::sin(modeA(1),modeK(1),j_);
                        valI = CartesianExactStateIds::sin(modeA(2),modeK(2),i_);

                     } else if(typeId == CartesianExactStateIds::COSSINCOS)
                     {
                        valK = CartesianExactStateIds::cos(modeA(0),modeK(0),k_);
                        valJ = CartesianExactStateIds::sin(modeA(1),modeK(1),j_);
                        valI = CartesianExactStateIds::cos(modeA(2),modeK(2),i_);

                     } else if(typeId == CartesianExactStateIds::SINCOSSIN)
                     {
                        valK = CartesianExactStateIds::sin(modeA(0),modeK(0),k_);
                        valJ = CartesianExactStateIds::cos(modeA(1),modeK(1),j_);
                        valI = CartesianExactStateIds::sin(modeA(2),modeK(2),i_);

                     } else if(typeId == CartesianExactStateIds::SINSINCOS)
                     {
                        valK = CartesianExactStateIds::sin(modeA(0),modeK(0),k_);
                        valJ = CartesianExactStateIds::sin(modeA(1),modeK(1),j_);
                        valI = CartesianExactStateIds::cos(modeA(2),modeK(2),i_);

                     } else if(typeId == CartesianExactStateIds::COSCOSSIN)
                     {
                        valK = CartesianExactStateIds::cos(modeA(0),modeK(0),k_);
                        valJ = CartesianExactStateIds::cos(modeA(1),modeK(1),j_);
                        valI = CartesianExactStateIds::sin(modeA(2),modeK(2),i_);

                     } else if(typeId == CartesianExactStateIds::SINCOSCOS)
                     {
                        valK = CartesianExactStateIds::sin(modeA(0),modeK(0),k_);
                        valJ = CartesianExactStateIds::cos(modeA(1),modeK(1),j_);
                        valI = CartesianExactStateIds::cos(modeA(2),modeK(2),i_);

                     } else if(typeId == CartesianExactStateIds::COSSINSIN)
                     {
                        valK = CartesianExactStateIds::cos(modeA(0),modeK(0),k_);
                        valJ = CartesianExactStateIds::sin(modeA(1),modeK(1),j_);
                        valI = CartesianExactStateIds::sin(modeA(2),modeK(2),i_);

                     } else
                     {
                        throw Exception("Unknown exact state");
                     }

                  } else if(typeId == CartesianExactStateIds::SPECIAL1)
                  {
                     valK = 3.31 + CartesianExactStateIds::poly(modeA(0),modeK(0),k_);
                     valJ = 2.71 + CartesianExactStateIds::cos(modeA(1),modeK(1),j_) + CartesianExactStateIds::sin(modeA(1),modeK(1),j_);
                     valI = -1.3 + CartesianExactStateIds::cos(modeA(2),modeK(2),i_) + CartesianExactStateIds::sin(modeA(2),modeK(2),i_);

                  } else if(typeId == CartesianExactStateIds::SPECIAL2)
                  {
                     valK = 3.31 + CartesianExactStateIds::poly(modeA(0),modeK(0),k_);
                     valJ = 2.71 + CartesianExactStateIds::cos(modeA(1),modeK(1),j_) - CartesianExactStateIds::sin(modeA(1),modeK(1),j_);
                     valI = -1.3 + CartesianExactStateIds::cos(modeA(2),modeK(2),i_) - CartesianExactStateIds::sin(modeA(2),modeK(2),i_);

                  } else if(typeId == CartesianExactStateIds::SPECIAL3)
                  {
                     valK = 3.31 + CartesianExactStateIds::poly(modeA(0),modeK(0),k_);
                     valJ = 2.71 + CartesianExactStateIds::cos(modeA(1),modeK(1),j_) + CartesianExactStateIds::sin(modeA(1),modeK(1),j_);
                     valI = -1.3 + CartesianExactStateIds::cos(modeA(2),modeK(2),i_) + CartesianExactStateIds::sin(modeA(2),modeK(2),i_);

                  } else
                  {
                     throw Exception("Unknown exact state");
                  }

                  rNLComp.setPoint(valK*valJ*valI, iI, iJ, iK);
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

}
}
