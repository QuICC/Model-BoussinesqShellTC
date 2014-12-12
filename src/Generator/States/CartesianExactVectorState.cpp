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

   void CartesianExactVectorState::setStateType(const FieldComponents::Physical::Id compId, const CartesianExactVectorState::StateTypeId id)
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
      StateTypeId typeId = this->mTypeId.find(compId)->second;
      Array modeA = this->mModeA.find(compId)->second;
      Array modeK = this->mModeK.find(compId)->second;

      if(typeId == CONSTANT)
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
                     if(typeId == POLYPOLYPOLY)
                     {
                        valK = this->poly(modeA(0),modeK(0),k_);
                        valJ = this->poly(modeA(1),modeK(1),j_);
                        valI = this->poly(modeA(2),modeK(2),i_);

                     } else
                     {
                        throw Exception("Unknown exact state");
                     }

                     // Generate solutions for TFT geometries
                  } else if(static_cast<int>(typeId) > 19 && static_cast<int>(typeId) < 30)
                  {
                     if(typeId == POLYCOSPOLY)
                     {
                        valK = this->poly(modeA(0),modeK(0),k_);
                        valJ = this->cos(modeA(1),modeK(1),j_);
                        valI = this->poly(modeA(2),modeK(2),i_);

                     } else if(typeId == POLYSINPOLY)
                     {
                        valK = this->poly(modeA(0),modeK(0),k_);
                        valJ = this->sin(modeA(1),modeK(1),j_);
                        valI = this->poly(modeA(2),modeK(2),i_);

                     } else
                     {
                        throw Exception("Unknown exact state");
                     }

                     // Generate solutions for TFF geometries
                  } else if(static_cast<int>(typeId) > 29 && static_cast<int>(typeId) < 50)
                  {
                     if(typeId == POLYCOSCOS)
                     {
                        valK = this->poly(modeA(0),modeK(0),k_);
                        valJ = this->cos(modeA(1),modeK(1),j_);
                        valI = this->cos(modeA(2),modeK(2),i_);

                     } else if(typeId == POLYSINSIN)
                     {
                        valK = this->poly(modeA(0),modeK(0),k_);
                        valJ = this->sin(modeA(1),modeK(1),j_);
                        valI = this->sin(modeA(2),modeK(2),i_);

                     } else if(typeId == POLYSINCOS)
                     {
                        valK = this->poly(modeA(0),modeK(0),k_);
                        valJ = this->sin(modeA(1),modeK(1),j_);
                        valI = this->cos(modeA(2),modeK(2),i_);

                     } else if(typeId == POLYCOSSIN)
                     {
                        valK = this->poly(modeA(0),modeK(0),k_);
                        valJ = this->cos(modeA(1),modeK(1),j_);
                        valI = this->sin(modeA(2),modeK(2),i_);

                     } else
                     {
                        throw Exception("Unknown exact state");
                     }

                     // Generate solutions for FFF geometries
                  } else if(static_cast<int>(typeId) > 39 && static_cast<int>(typeId) < 50)
                  {
                     if(typeId == COSCOSCOS)
                     {
                        valK = this->cos(modeA(0),modeK(0),k_);
                        valJ = this->cos(modeA(1),modeK(1),j_);
                        valI = this->cos(modeA(2),modeK(2),i_);

                     } else if(typeId == SINSINSIN)
                     {
                        valK = this->sin(modeA(0),modeK(0),k_);
                        valJ = this->sin(modeA(1),modeK(1),j_);
                        valI = this->sin(modeA(2),modeK(2),i_);

                     } else if(typeId == COSSINCOS)
                     {
                        valK = this->cos(modeA(0),modeK(0),k_);
                        valJ = this->sin(modeA(1),modeK(1),j_);
                        valI = this->cos(modeA(2),modeK(2),i_);

                     } else if(typeId == SINCOSSIN)
                     {
                        valK = this->sin(modeA(0),modeK(0),k_);
                        valJ = this->cos(modeA(1),modeK(1),j_);
                        valI = this->sin(modeA(2),modeK(2),i_);

                     } else if(typeId == SINSINCOS)
                     {
                        valK = this->sin(modeA(0),modeK(0),k_);
                        valJ = this->sin(modeA(1),modeK(1),j_);
                        valI = this->cos(modeA(2),modeK(2),i_);

                     } else if(typeId == COSCOSSIN)
                     {
                        valK = this->cos(modeA(0),modeK(0),k_);
                        valJ = this->cos(modeA(1),modeK(1),j_);
                        valI = this->sin(modeA(2),modeK(2),i_);

                     } else if(typeId == SINCOSCOS)
                     {
                        valK = this->sin(modeA(0),modeK(0),k_);
                        valJ = this->cos(modeA(1),modeK(1),j_);
                        valI = this->cos(modeA(2),modeK(2),i_);

                     } else if(typeId == COSSINSIN)
                     {
                        valK = this->cos(modeA(0),modeK(0),k_);
                        valJ = this->sin(modeA(1),modeK(1),j_);
                        valI = this->sin(modeA(2),modeK(2),i_);

                     } else
                     {
                        throw Exception("Unknown exact state");
                     }

                  } else if(typeId == SPECIAL1)
                  {
                     valK = 3.31 + this->poly(modeA(0),modeK(0),k_);
                     valJ = 2.71 + this->cos(modeA(1),modeK(1),j_) + this->sin(modeA(1),modeK(1),j_);
                     valI = -1.3 + this->cos(modeA(2),modeK(2),i_) + this->sin(modeA(2),modeK(2),i_);

                  } else if(typeId == SPECIAL2)
                  {
                     valK = 3.31 + this->poly(modeA(0),modeK(0),k_);
                     valJ = 2.71 + this->cos(modeA(1),modeK(1),j_) - this->sin(modeA(1),modeK(1),j_);
                     valI = -1.3 + this->cos(modeA(2),modeK(2),i_) - this->sin(modeA(2),modeK(2),i_);

                  } else if(typeId == SPECIAL3)
                  {
                     valK = 3.31 + this->poly(modeA(0),modeK(0),k_);
                     valJ = 2.71 + this->cos(modeA(1),modeK(1),j_) + this->sin(modeA(1),modeK(1),j_);
                     valI = -1.3 + this->cos(modeA(2),modeK(2),i_) + this->sin(modeA(2),modeK(2),i_);

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

   MHDFloat CartesianExactVectorState::cos(const MHDFloat amplitude, const MHDFloat mode, const MHDFloat theta) const
   {
      return amplitude*std::cos(mode*theta);
   }

   MHDFloat CartesianExactVectorState::sin(const MHDFloat amplitude, const MHDFloat mode, const MHDFloat theta) const
   {
      return amplitude*std::sin(mode*theta);
   }

   MHDFloat CartesianExactVectorState::poly(const MHDFloat amplitude, const MHDFloat mode, const MHDFloat x) const
   {
      MHDFloat val;

      if(mode == CartesianExactVectorState::PCOS)
      {
         val = this->cos(amplitude,mode,Math::PI*(x-1)/2.0);
      } else if(mode == CartesianExactVectorState::PSIN)
      {
         val = this->sin(amplitude,mode,Math::PI*(x-1)/2.0);
      } else
      {
         val = amplitude*std::pow(x,mode);
      }

      return val;
   }

}
}
