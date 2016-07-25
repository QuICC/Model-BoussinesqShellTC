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

   void CartesianExactScalarState::setModeOptions(const MHDFloat a1, const MHDFloat k1, const MHDFloat a2, const MHDFloat k2)
   {
      this->mModeA.resize(2);
      this->mModeA(0) = a1;
      this->mModeA(1) = a2;

      this->mModeK.resize(2);
      this->mModeK(0) = k1;
      this->mModeK(1) = k2;
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

      if(this->mTypeId == CartesianExactStateIds::CONSTANT)
      {
         rNLComp.rData().setConstant(this->mModeA.prod());

      } else if(this->mTypeId == CartesianExactStateIds::PEYRET1DA)
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

                  MHDFloat val = (-1.0 + 2.0*Math::PI)*std::cos(Math::PI*k_/2.0)/(2.0*Math::PI);

                  rNLComp.setPoint(val, iI, iJ, iK);
               }
            }
         }
      } else
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
                  if(static_cast<int>(this->mTypeId) < 100)
                  {
                     Array grid(3);
                     grid(0) = k_;
                     grid(1) = j_;
                     grid(2) = i_;
                     val = CartesianExactStateIds::exact3D(this->mTypeId, this->mModeA, this->mModeK, grid);
                  } else if(static_cast<int>(this->mTypeId) >= 100)
                  {
                     Array grid(2);
                     grid(0) = j_;
                     grid(1) = i_;
                     val = CartesianExactStateIds::exact2D(this->mTypeId, this->mModeA, this->mModeK, grid);
                  }

                  rNLComp.setPoint(val, iI, iJ, iK);
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

      // Add unknown to requirements: is scalar?, need spectral?, need physical?, need diff?, need curl?, need diff2?
      this->mRequirements.addField(this->name(), FieldRequirement(true, true, true, false));
   }

}
}
