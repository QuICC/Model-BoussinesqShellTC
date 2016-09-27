/** 
 * @file CylinderExactVectorState.cpp
 * @brief Source of the implementation of the equation to generate an exact vector solution in a cylinder
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
#include "Generator/States/CylinderExactVectorState.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "TypeSelectors/TransformSelector.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   CylinderExactVectorState::CylinderExactVectorState(SharedEquationParameters spEqParams)
      : IVectorEquation(spEqParams)
   {
   }

   CylinderExactVectorState::~CylinderExactVectorState()
   {
   }

   void CylinderExactVectorState::setIdentity(const PhysicalNames::Id name)
   {
      // Set the name
      this->setName(name);

      // Set the variable requirements
      this->setRequirements();
   }

   void CylinderExactVectorState::setStateType(const FieldComponents::Physical::Id compId, const CylinderExactStateIds::Id id)
   {
      this->mTypeId.insert(std::make_pair(compId, id));
   }

   void CylinderExactVectorState::setModeOptions(const FieldComponents::Physical::Id compId, const MHDFloat a1, const MHDFloat k1, const MHDFloat a2, const MHDFloat k2, const MHDFloat a3, const MHDFloat k3)
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

   void CylinderExactVectorState::setCoupling()
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

   void CylinderExactVectorState::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id compId) const
   {
      CylinderExactStateIds::Id typeId = this->mTypeId.find(compId)->second;
      Array modeA = this->mModeA.find(compId)->second;
      Array modeK = this->mModeK.find(compId)->second;

      if(typeId == CylinderExactStateIds::CONSTANT)
      {
         rNLComp.rData().setConstant(modeA(0)*modeA(1)*modeA(2));
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

                  if(typeId == CylinderExactStateIds::POLYCOSPOLY)
                  {
                     MHDFloat valR = CylinderExactStateIds::poly(modeA(0),modeK(0),r_);
                     MHDFloat valT = CylinderExactStateIds::cos(modeA(1),modeK(1),t_);
                     MHDFloat valZ = CylinderExactStateIds::poly(modeA(2),modeK(2),z_);

                     rNLComp.setPoint(valR*valT*valZ, iZ, iT, iR);

                  } else if(typeId == CylinderExactStateIds::POLYSINPOLY)
                  {
                     MHDFloat valR = CylinderExactStateIds::poly(modeA(0),modeK(0),r_);
                     MHDFloat valT = CylinderExactStateIds::sin(modeA(1),modeK(1),t_);
                     MHDFloat valZ = CylinderExactStateIds::poly(modeA(2),modeK(2),z_);

                     rNLComp.setPoint(valR*valT*valZ, iZ, iT, iR);

                  } else if(typeId == CylinderExactStateIds::TESTUNITSPECTRUM)
                  {
                     MHDFloat val = 0.0;
                     MHDFloat valR = 0.0;
                     MHDFloat valT = 0.0;
                     MHDFloat valZ = 0.0;
                     for(int jM = 1; jM <= modeK(1);  ++jM)
                     {
                        valR = 0.0;
                        for(int jR = 0; jR <= modeK(0); jR++)
                        {
                           valR += std::pow(r_, 2*jR);
                        }
                        valR *= std::pow(r_, jM);

                        valT = CylinderExactStateIds::cos(1.0,jM,t_) + CylinderExactStateIds::sin(1.0,jM,t_);

                        valZ = 0.0;
                        for(int jZ = 0; jZ <= modeK(2); jZ++)
                        {
                           valZ += CylinderExactStateIds::chebyshev(1.0, jZ, z_);
                        }
                        val += valR*valT*valZ;
                     }

                     rNLComp.setPoint(val, iZ, iT, iR);

                  } else
                  {
                     throw Exception("Unknown exact state");
                  }
               }
            }
         }
      }
   }

    Datatypes::SpectralScalarType::PointType CylinderExactVectorState::sourceTerm(FieldComponents::Spectral::Id compId, const int i, const int j, const int k) const
    {
      return Datatypes::SpectralScalarType::PointType(0);
    }

   void CylinderExactVectorState::setRequirements()
   {
      // Set solver timing
      this->setSolveTiming(SolveTiming::AFTER);

      // Add unknown to requirements: is scalar?, need spectral?, need physical?, need diff?(, need curl?)
      this->mRequirements.addField(this->name(), FieldRequirement(false, true, true, false));
   }

   void CylinderExactVectorState::setNLComponents()
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
