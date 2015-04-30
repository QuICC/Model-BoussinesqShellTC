/** 
 * @file AnnulusExactScalarState.cpp
 * @brief Source of the implementation of the equation to generate an exact scalar solution in an annulus
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
#include "Generator/States/AnnulusExactScalarState.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "TypeSelectors/TransformSelector.hpp"
#include "TypeSelectors/EquationEigenSelector.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   AnnulusExactScalarState::AnnulusExactScalarState(SharedEquationParameters spEqParams)
      : IScalarEquation(spEqParams), mTypeId(AnnulusExactStateIds::CONSTANT), mModeA(3), mModeK(3)
   {
   }

   AnnulusExactScalarState::~AnnulusExactScalarState()
   {
   }

   void AnnulusExactScalarState::setIdentity(const PhysicalNames::Id name)
   {
      // Set the name
      this->setName(name);

      // Set the variable requirements
      this->setRequirements();
   }

   void AnnulusExactScalarState::setStateType(const AnnulusExactStateIds::Id id)
   {
      this->mTypeId = id;
   }

   void AnnulusExactScalarState::setModeOptions(const MHDFloat a1, const MHDFloat k1, const MHDFloat a2, const MHDFloat k2, const MHDFloat a3, const MHDFloat k3)
   {
      this->mModeA(0) = a1;
      this->mModeA(1) = a2;
      this->mModeA(2) = a3;

      this->mModeK(0) = k1;
      this->mModeK(1) = k2;
      this->mModeK(2) = k3;
   }

   void AnnulusExactScalarState::setCoupling()
   {
      this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::TRIVIAL, 0, true, false, false);
   }

   void AnnulusExactScalarState::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id compId) const
   {
      // Assert on scalar component is used
      assert(compId == FieldComponents::Physical::SCALAR);

      if(this->mTypeId == AnnulusExactStateIds::CONSTANT)
      {
         rNLComp.rData().setConstant(this->mModeA(0)*this->mModeA(1)*this->mModeA(2));
      } else
      {
         int nR = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D,Dimensions::Space::PHYSICAL);
         int nT = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM2D,Dimensions::Space::PHYSICAL);
         int nZ = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM3D,Dimensions::Space::PHYSICAL);

         Array gR = Transform::TransformSelector<Dimensions::Transform::TRA1D>::Type::generateGrid(nR, 1, 0.35);
         Array gT = Transform::TransformSelector<Dimensions::Transform::TRA2D>::Type::generateGrid(nT);
         Array gZ = Transform::TransformSelector<Dimensions::Transform::TRA3D>::Type::generateGrid(nZ);

         MHDFloat r_;
         MHDFloat t_;
         MHDFloat z_;

         nR = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();
         for(int iR = 0; iR < nR; ++iR)
         {
            r_ = gR(this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(iR));
            nT = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(iR);
            for(int iT = 0; iT < nT; ++iT)
            {
               t_ = gT(this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT2D>(iT, iR));
               for(int iZ = 0; iZ < nZ; ++iZ)
               {
                  z_ = gZ(iZ);

                  MHDFloat valZ = 0.0;
                  MHDFloat valT = 0.0;
                  MHDFloat valR = 0.0;

                  if(this->mTypeId == AnnulusExactStateIds::POLYCOSPOLY)
                  {
                     valR = AnnulusExactStateIds::poly(this->mModeA(0),this->mModeK(0),r_);
                     valT = AnnulusExactStateIds::cos(this->mModeA(1),this->mModeK(1),t_);
                     valZ = AnnulusExactStateIds::poly(this->mModeA(2),this->mModeK(2),z_);

                  } else if(this->mTypeId == AnnulusExactStateIds::POLYSINPOLY)
                  {
                     valR = AnnulusExactStateIds::poly(this->mModeA(0),this->mModeK(0),r_);
                     valT = AnnulusExactStateIds::sin(this->mModeA(1),this->mModeK(1),t_);
                     valZ = AnnulusExactStateIds::poly(this->mModeA(2),this->mModeK(2),z_);

                  } else
                  {
                     throw Exception("Unknown exact state");
                  }

                  rNLComp.setPoint(valR*valT*valZ, iZ, iT, iR);
               }
            }
         }
      }
   }

    Datatypes::SpectralScalarType::PointType AnnulusExactScalarState::sourceTerm(FieldComponents::Spectral::Id compId, const int i, const int j, const int k) const
    {
      // Assert on scalar component is used
      assert(compId == FieldComponents::Spectral::SCALAR);

      return Datatypes::SpectralScalarType::PointType(0);
    }

   void AnnulusExactScalarState::setRequirements()
   {
      // Set solver timing
      this->setSolveTiming(SolveTiming::AFTER);

      // Add unknown to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(this->name(), FieldRequirement(true, true, false, false));
   }

}
}
