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
#include "TypeSelectors/TransformSelector.hpp"
#include "TypeSelectors/EquationEigenSelector.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   CartesianExactScalarState::CartesianExactScalarState(SharedEquationParameters spEqParams)
      : IScalarEquation(spEqParams), mTypeId(CartesianExactStateIds::CONSTANT), mModeA(3), mModeK(3)
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

   void CartesianExactScalarState::setStateType(const CartesianExactStateIds::Id id)
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
      this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::TRIVIAL, 0, true, false, false, false);
   }

   void CartesianExactScalarState::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id compId) const
   {
      // Assert on scalar component is used
      assert(compId == FieldComponents::Physical::SCALAR);

      if(this->mTypeId == CartesianExactStateIds::CONSTANT)
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
                     if(this->mTypeId == CartesianExactStateIds::CartesianExactStateIds::POLYPOLYPOLY)
                     {
                        valK = CartesianExactStateIds::poly(this->mModeA(0),this->mModeK(0),k_);
                        valJ = CartesianExactStateIds::poly(this->mModeA(1),this->mModeK(1),j_);
                        valI = CartesianExactStateIds::poly(this->mModeA(2),this->mModeK(2),i_);

                     } else
                     {
                        throw Exception("Unknown exact state");
                     }

                     // Generate solutions for TFT geometries
                  } else if(static_cast<int>(this->mTypeId) > 19 && static_cast<int>(this->mTypeId) < 30)
                  {
                     if(this->mTypeId == CartesianExactStateIds::POLYCOSPOLY)
                     {
                        valK = CartesianExactStateIds::poly(this->mModeA(0),this->mModeK(0),k_);
                        valJ = CartesianExactStateIds::cos(this->mModeA(1),this->mModeK(1),j_);
                        valI = CartesianExactStateIds::poly(this->mModeA(2),this->mModeK(2),i_);

                     } else if(this->mTypeId == CartesianExactStateIds::POLYSINPOLY)
                     {
                        valK = CartesianExactStateIds::poly(this->mModeA(0),this->mModeK(0),k_);
                        valJ = CartesianExactStateIds::sin(this->mModeA(1),this->mModeK(1),j_);
                        valI = CartesianExactStateIds::poly(this->mModeA(2),this->mModeK(2),i_);

                     } else
                     {
                        throw Exception("Unknown exact state");
                     }

                     // Generate solutions for TFF geometries
                  } else if(static_cast<int>(this->mTypeId) > 29 && static_cast<int>(this->mTypeId) < 50)
                  {
                     if(this->mTypeId == CartesianExactStateIds::POLYCOSCOS)
                     {
                        valK = CartesianExactStateIds::poly(this->mModeA(0),this->mModeK(0),k_);
                        valJ = CartesianExactStateIds::cos(this->mModeA(1),this->mModeK(1),j_);
                        valI = CartesianExactStateIds::cos(this->mModeA(2),this->mModeK(2),i_);

                     } else if(this->mTypeId == CartesianExactStateIds::POLYSINSIN)
                     {
                        valK = CartesianExactStateIds::poly(this->mModeA(0),this->mModeK(0),k_);
                        valJ = CartesianExactStateIds::sin(this->mModeA(1),this->mModeK(1),j_);
                        valI = CartesianExactStateIds::sin(this->mModeA(2),this->mModeK(2),i_);

                     } else if(this->mTypeId == CartesianExactStateIds::POLYSINCOS)
                     {
                        valK = CartesianExactStateIds::poly(this->mModeA(0),this->mModeK(0),k_);
                        valJ = CartesianExactStateIds::sin(this->mModeA(1),this->mModeK(1),j_);
                        valI = CartesianExactStateIds::cos(this->mModeA(2),this->mModeK(2),i_);

                     } else if(this->mTypeId == CartesianExactStateIds::POLYCOSSIN)
                     {
                        valK = CartesianExactStateIds::poly(this->mModeA(0),this->mModeK(0),k_);
                        valJ = CartesianExactStateIds::cos(this->mModeA(1),this->mModeK(1),j_);
                        valI = CartesianExactStateIds::sin(this->mModeA(2),this->mModeK(2),i_);

                     } else
                     {
                        throw Exception("Unknown exact state");
                     }

                     // Generate solutions for FFF geometries
                  } else if(static_cast<int>(this->mTypeId) > 39 && static_cast<int>(this->mTypeId) < 50)
                  {
                     if(this->mTypeId == CartesianExactStateIds::COSCOSCOS)
                     {
                        valK = CartesianExactStateIds::cos(this->mModeA(0),this->mModeK(0),k_);
                        valJ = CartesianExactStateIds::cos(this->mModeA(1),this->mModeK(1),j_);
                        valI = CartesianExactStateIds::cos(this->mModeA(2),this->mModeK(2),i_);

                     } else if(this->mTypeId == CartesianExactStateIds::SINSINSIN)
                     {
                        valK = CartesianExactStateIds::sin(this->mModeA(0),this->mModeK(0),k_);
                        valJ = CartesianExactStateIds::sin(this->mModeA(1),this->mModeK(1),j_);
                        valI = CartesianExactStateIds::sin(this->mModeA(2),this->mModeK(2),i_);

                     } else if(this->mTypeId == CartesianExactStateIds::COSSINCOS)
                     {
                        valK = CartesianExactStateIds::cos(this->mModeA(0),this->mModeK(0),k_);
                        valJ = CartesianExactStateIds::sin(this->mModeA(1),this->mModeK(1),j_);
                        valI = CartesianExactStateIds::cos(this->mModeA(2),this->mModeK(2),i_);

                     } else if(this->mTypeId == CartesianExactStateIds::SINCOSSIN)
                     {
                        valK = CartesianExactStateIds::sin(this->mModeA(0),this->mModeK(0),k_);
                        valJ = CartesianExactStateIds::cos(this->mModeA(1),this->mModeK(1),j_);
                        valI = CartesianExactStateIds::sin(this->mModeA(2),this->mModeK(2),i_);

                     } else if(this->mTypeId == CartesianExactStateIds::SINSINCOS)
                     {
                        valK = CartesianExactStateIds::sin(this->mModeA(0),this->mModeK(0),k_);
                        valJ = CartesianExactStateIds::sin(this->mModeA(1),this->mModeK(1),j_);
                        valI = CartesianExactStateIds::cos(this->mModeA(2),this->mModeK(2),i_);

                     } else if(this->mTypeId == CartesianExactStateIds::COSCOSSIN)
                     {
                        valK = CartesianExactStateIds::cos(this->mModeA(0),this->mModeK(0),k_);
                        valJ = CartesianExactStateIds::cos(this->mModeA(1),this->mModeK(1),j_);
                        valI = CartesianExactStateIds::sin(this->mModeA(2),this->mModeK(2),i_);

                     } else if(this->mTypeId == CartesianExactStateIds::SINCOSCOS)
                     {
                        valK = CartesianExactStateIds::sin(this->mModeA(0),this->mModeK(0),k_);
                        valJ = CartesianExactStateIds::cos(this->mModeA(1),this->mModeK(1),j_);
                        valI = CartesianExactStateIds::cos(this->mModeA(2),this->mModeK(2),i_);

                     } else if(this->mTypeId == CartesianExactStateIds::COSSINSIN)
                     {
                        valK = CartesianExactStateIds::cos(this->mModeA(0),this->mModeK(0),k_);
                        valJ = CartesianExactStateIds::sin(this->mModeA(1),this->mModeK(1),j_);
                        valI = CartesianExactStateIds::sin(this->mModeA(2),this->mModeK(2),i_);

                     } else
                     {
                        throw Exception("Unknown exact state");
                     }

                  } else if(this->mTypeId == CartesianExactStateIds::SPECIAL1)
                  {
                     valK = 3.31 + CartesianExactStateIds::poly(this->mModeA(0),this->mModeK(0),k_);
                     valJ = 2.71 + CartesianExactStateIds::cos(this->mModeA(1),this->mModeK(1),j_) + CartesianExactStateIds::sin(this->mModeA(1),this->mModeK(1),j_);
                     valI = -1.3 + CartesianExactStateIds::cos(this->mModeA(2),this->mModeK(2),i_) + CartesianExactStateIds::sin(this->mModeA(2),this->mModeK(2),i_);

                  } else if(this->mTypeId == CartesianExactStateIds::SPECIAL2)
                  {
                     valK = 3.31 + CartesianExactStateIds::poly(this->mModeA(0),this->mModeK(0),k_);
                     valJ = 2.71 + CartesianExactStateIds::cos(this->mModeA(1),this->mModeK(1),j_) - CartesianExactStateIds::sin(this->mModeA(1),this->mModeK(1),j_);
                     valI = -1.3 + CartesianExactStateIds::cos(this->mModeA(2),this->mModeK(2),i_) - CartesianExactStateIds::sin(this->mModeA(2),this->mModeK(2),i_);

                  } else if(this->mTypeId == CartesianExactStateIds::SPECIAL3)
                  {
                     valK = 3.31 + CartesianExactStateIds::poly(this->mModeA(0),this->mModeK(0),k_);
                     valJ = 2.71 + CartesianExactStateIds::cos(this->mModeA(1),this->mModeK(1),j_) + CartesianExactStateIds::sin(this->mModeA(1),this->mModeK(1),j_);
                     valI = -1.3 + CartesianExactStateIds::cos(this->mModeA(2),this->mModeK(2),i_) + CartesianExactStateIds::sin(this->mModeA(2),this->mModeK(2),i_);

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
      // Set solver timing
      this->setSolveTiming(SolveTiming::AFTER);

      // Add unknown to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(this->name(), FieldRequirement(true, true, true, false));
   }

}
}
