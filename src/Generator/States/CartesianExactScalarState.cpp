/** 
 * @file CartesianExactScalarState.cpp
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
#include "Generator/States/CartesianExactScalarState.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "TypeSelectors/SpectralOperatorSelector.hpp"
#include "TypeSelectors/TransformSelector.hpp"
#include "TypeSelectors/EquationEigenSelector.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   CartesianExactScalarState::CartesianExactScalarState(SharedEquationParameters spEqParams)
      : IScalarEquation(spEqParams), mTypeId(CONSTANT), mModeA(3), mModeK(3)
   {
   }

   CartesianExactScalarState::~CartesianExactScalarState()
   {
   }

   void CartesianExactScalarState::setIdentity(const PhysicalNames::Id name)
   {
      // Set the name
      this->setName(name);

      // Set the variable requirements
      this->setRequirements();
   }

   void CartesianExactScalarState::setStateType(const CartesianExactScalarState::StateTypeId id)
   {
      this->mTypeId = id;
   }

   void CartesianExactScalarState::setModeOptions(const MHDFloat a1, const MHDFloat k1, const MHDFloat a2, const MHDFloat k2, const MHDFloat a3, const MHDFloat k3)
   {
      this->mModeA(0) = a1;
      this->mModeA(1) = a2;
      this->mModeA(2) = a3;

      this->mModeK(0) = k1;
      this->mModeK(1) = k2;
      this->mModeK(2) = k3;
   }

   void CartesianExactScalarState::setCoupling()
   {
      this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::TRIVIAL, 0, true, false, false);
   }

   void CartesianExactScalarState::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id compId) const
   {
      // Assert on scalar component is used
      assert(compId == FieldComponents::Physical::SCALAR);

      if(this->mTypeId == CONSTANT)
      {
         rNLComp.rData().setConstant(this->mModeA(0)*this->mModeA(1)*this->mModeA(2));
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
                  if(static_cast<int>(this->mTypeId) > 9 && static_cast<int>(this->mTypeId) < 20)
                  {
                     if(this->mTypeId == POLYPOLYPOLY)
                     {
                        valK = this->poly(0,k_);
                        valJ = this->poly(1,j_);
                        valI = this->poly(2,i_);

                     } else
                     {
                        throw Exception("Unknown exact state");
                     }

                     // Generate solutions for TFT geometries
                  } else if(static_cast<int>(this->mTypeId) > 19 && static_cast<int>(this->mTypeId) < 30)
                  {
                     if(this->mTypeId == POLYCOSPOLY)
                     {
                        valK = this->poly(0,k_);
                        valJ = this->cos(1,j_);
                        valI = this->poly(2,i_);

                     } else if(this->mTypeId == POLYSINPOLY)
                     {
                        valK = this->poly(0,k_);
                        valJ = this->sin(1,j_);
                        valI = this->poly(2,i_);

                     } else
                     {
                        throw Exception("Unknown exact state");
                     }

                     // Generate solutions for TFF geometries
                  } else if(static_cast<int>(this->mTypeId) > 29 && static_cast<int>(this->mTypeId) < 40)
                  {
                     if(this->mTypeId == POLYCOSCOS)
                     {
                        valK = this->poly(0,k_);
                        valJ = this->cos(1,j_);
                        valI = this->cos(2,i_);

                     } else if(this->mTypeId == POLYSINSIN)
                     {
                        valK = this->poly(0,k_);
                        valJ = this->sin(1,j_);
                        valI = this->sin(2,i_);

                     } else if(this->mTypeId == POLYSINCOS)
                     {
                        valK = this->poly(0,k_);
                        valJ = this->sin(1,j_);
                        valI = this->cos(2,i_);

                     } else if(this->mTypeId == POLYCOSSIN)
                     {
                        valK = this->poly(0,k_);
                        valJ = this->cos(1,j_);
                        valI = this->sin(2,i_);

                     } else if(this->mTypeId == PSINCOSSIN)
                     {
                        valK = this->sin(0,Math::PI*(k_-1)/2.0);
                        valJ = this->cos(1,j_);
                        valI = this->sin(2,i_);

                     } else if(this->mTypeId == PCOSCOSSIN)
                     {
                        valK = this->cos(0,Math::PI*(k_-1)/2.0);
                        valJ = this->cos(1,j_);
                        valI = this->sin(2,i_);
                     } else
                     {
                        throw Exception("Unknown exact state");
                     }

                     // Generate solutions for FFF geometries
                  } else if(static_cast<int>(this->mTypeId) > 39 && static_cast<int>(this->mTypeId) < 50)
                  {
                     if(this->mTypeId == COSCOSCOS)
                     {
                        valK = this->cos(0,k_);
                        valJ = this->cos(1,j_);
                        valI = this->cos(2,i_);

                     } else if(this->mTypeId == SINSINSIN)
                     {
                        valK = this->sin(0,k_);
                        valJ = this->sin(1,j_);
                        valI = this->sin(2,i_);

                     } else if(this->mTypeId == COSSINCOS)
                     {
                        valK = this->cos(0,k_);
                        valJ = this->sin(1,j_);
                        valI = this->cos(2,i_);

                     } else if(this->mTypeId == SINCOSSIN)
                     {
                        valK = this->sin(0,k_);
                        valJ = this->cos(1,j_);
                        valI = this->sin(2,i_);

                     } else if(this->mTypeId == SINSINCOS)
                     {
                        valK = this->sin(0,k_);
                        valJ = this->sin(1,j_);
                        valI = this->cos(2,i_);

                     } else if(this->mTypeId == COSCOSSIN)
                     {
                        valK = this->cos(0,k_);
                        valJ = this->cos(1,j_);
                        valI = this->sin(2,i_);

                     } else if(this->mTypeId == SINCOSCOS)
                     {
                        valK = this->sin(0,k_);
                        valJ = this->cos(1,j_);
                        valI = this->cos(2,i_);

                     } else if(this->mTypeId == COSSINSIN)
                     {
                        valK = this->cos(0,k_);
                        valJ = this->sin(1,j_);
                        valI = this->sin(2,i_);

                     } else
                     {
                        throw Exception("Unknown exact state");
                     }
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

    Datatypes::SpectralScalarType::PointType CartesianExactScalarState::sourceTerm(FieldComponents::Spectral::Id compId, const int i, const int j, const int k) const
    {
      // Assert on scalar component is used
      assert(compId == FieldComponents::Spectral::SCALAR);

      return Datatypes::SpectralScalarType::PointType(0);
    }

   void CartesianExactScalarState::setRequirements()
   {
      // Add unknown to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(this->name(), FieldRequirement(true, true, true, false));
   }

   MHDFloat CartesianExactScalarState::cos(const int idx, const MHDFloat theta) const
   {
      return this->mModeA(idx)*std::cos(this->mModeK(idx)*theta);
   }

   MHDFloat CartesianExactScalarState::sin(const int idx, const MHDFloat theta) const
   {
      return this->mModeA(idx)*std::sin(this->mModeK(idx)*theta);
   }

   MHDFloat CartesianExactScalarState::poly(const int idx, const MHDFloat x) const
   {
      return this->mModeA(idx)*std::pow(x,this->mModeK(idx));
   }

}
}
