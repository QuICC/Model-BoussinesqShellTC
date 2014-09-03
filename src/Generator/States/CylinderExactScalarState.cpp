/** 
 * @file CylinderExactScalarState.cpp
 * @brief Source of the implementation of the equation to generate an exact scalar solution in a cylinder
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
#include "Generator/States/CylinderExactScalarState.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "TypeSelectors/TransformSelector.hpp"
#include "TypeSelectors/EquationEigenSelector.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   CylinderExactScalarState::CylinderExactScalarState(SharedEquationParameters spEqParams)
      : IScalarEquation(spEqParams), mTypeId(CONSTANT), mModeA(3), mModeK(3)
   {
   }

   CylinderExactScalarState::~CylinderExactScalarState()
   {
   }

   void CylinderExactScalarState::setIdentity(const PhysicalNames::Id name)
   {
      // Set the name
      this->setName(name);

      // Set the variable requirements
      this->setRequirements();
   }

   void CylinderExactScalarState::setStateType(const CylinderExactScalarState::StateTypeId id)
   {
      this->mTypeId = id;
   }

   void CylinderExactScalarState::setModeOptions(const MHDFloat a1, const MHDFloat k1, const MHDFloat a2, const MHDFloat k2, const MHDFloat a3, const MHDFloat k3)
   {
      this->mModeA(0) = a1;
      this->mModeA(1) = a2;
      this->mModeA(2) = a3;

      this->mModeK(0) = k1;
      this->mModeK(1) = k2;
      this->mModeK(2) = k3;
   }

   void CylinderExactScalarState::setCoupling()
   {
      this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::TRIVIAL, 0, true, false, false, false);
   }

   void CylinderExactScalarState::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id compId) const
   {
      // Assert on scalar component is used
      assert(compId == FieldComponents::Physical::SCALAR);

      if(this->mTypeId == CONSTANT)
      {
         rNLComp.rData().setConstant(this->mModeA(0)*this->mModeA(1)*this->mModeA(2));
      } else
      {
         int nR = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D,Dimensions::Space::PHYSICAL);
         int nT = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM2D,Dimensions::Space::PHYSICAL);
         int nZ = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM3D,Dimensions::Space::PHYSICAL);

         Array gR = Transform::TransformSelector<Dimensions::Transform::TRA1D>::Type::generateGrid(nR);
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

                  if(this->mTypeId == POLYCOSPOLY)
                  {
                     valR = this->poly(0,r_);
                     valT = this->cos(1,t_);
                     valZ = this->poly(2,z_);

                  } else if(this->mTypeId == POLYSINPOLY)
                  {
                     valR = this->poly(0,r_);
                     valT = this->sin(1,t_);
                     valZ = this->poly(2,z_);

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

    Datatypes::SpectralScalarType::PointType CylinderExactScalarState::sourceTerm(FieldComponents::Spectral::Id compId, const int i, const int j, const int k) const
    {
      // Assert on scalar component is used
      assert(compId == FieldComponents::Spectral::SCALAR);

      return Datatypes::SpectralScalarType::PointType(0);
    }

   void CylinderExactScalarState::setRequirements()
   {
      // Set solver timing
      this->setSolveTiming(SolveTiming::AFTER);

      // Add unknown to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(this->name(), FieldRequirement(true, true, false, false));
   }

   MHDFloat CylinderExactScalarState::cos(const int idx, const MHDFloat theta) const
   {
      return this->mModeA(idx)*std::cos(this->mModeK(idx)*theta);
   }

   MHDFloat CylinderExactScalarState::sin(const int idx, const MHDFloat theta) const
   {
      return this->mModeA(idx)*std::sin(this->mModeK(idx)*theta);
   }

   MHDFloat CylinderExactScalarState::poly(const int idx, const MHDFloat x) const
   {
      MHDFloat val;

      if(this->mModeK(idx) == CylinderExactScalarState::PCOS)
      {
         val = this->cos(idx,Math::PI*(x-1)/2.0);
      } else if(this->mModeK(idx) == CylinderExactScalarState::PSIN)
      {
         val = this->sin(idx,Math::PI*(x-1)/2.0);
      } else
      {
         val = this->mModeA(idx)*std::pow(x,this->mModeK(idx));
      }

      return val;
   }

}
}
