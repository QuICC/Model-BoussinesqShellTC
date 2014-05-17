/** 
 * @file ExactScalarState.cpp
 * @brief Source of the implementation of the equation to generate an exact scalar solution
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
#include "Generator/States/ExactScalarState.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "TypeSelectors/SpectralOperatorSelector.hpp"
#include "TypeSelectors/TransformSelector.hpp"
#include "TypeSelectors/EquationEigenSelector.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   ExactScalarState::ExactScalarState(const std::string& pyName, SharedEquationParameters spEqParams)
      : IScalarEquation(pyName, spEqParams), mTypeId(CONSTANT), mSineA(2), mSineN(2)
   {
   }

   ExactScalarState::~ExactScalarState()
   {
   }

   void ExactScalarState::setIdentity(const PhysicalNames::Id name)
   {
      // Set the name
      this->setName(name);

      // Set the variable requirements
      this->setRequirements();
   }

   void ExactScalarState::setStateType(const ExactScalarState::StateTypeId id)
   {
      this->mTypeId = id;
   }

   void ExactScalarState::setSineOptions(const MHDFloat aX, const MHDFloat kX, const MHDFloat aZ, const MHDFloat kZ)
   {
      this->mSineA(0) = aX;
      this->mSineA(1) = aZ;

      this->mSineN(0) = kX;
      this->mSineN(1) = kZ;
   }

   void ExactScalarState::setCoupling()
   {
      this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::TRIVIAL, 0, true, false, false);
   }

   void ExactScalarState::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id compId) const
   {
      // Assert on scalar component is used
      assert(compId == FieldComponents::Physical::SCALAR);

      if(this->mTypeId == CONSTANT)
      {
         rNLComp.rData().setConstant(42.0);
      } else if(this->mTypeId == SINESINE)
      {
         int nX = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();
         int nY = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM2D,Dimensions::Space::PHYSICAL);
         int nZ = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM3D,Dimensions::Space::PHYSICAL);

         Array xGrid = Transform::TransformSelector<Dimensions::Transform::TRA1D>::Type::generateGrid(nX);
         Array yGrid = Transform::TransformSelector<Dimensions::Transform::TRA2D>::Type::generateGrid(nY);
         Array zGrid = Transform::TransformSelector<Dimensions::Transform::TRA3D>::Type::generateGrid(nZ);

         for(int iX = 0; iX < nX; ++iX)
         {
            for(int iY = 0; iY < nY; ++iY)
            {
               Array sinZ = this->mSineA(1)*(this->mSineN(1)*(Math::PI/2)*(1+zGrid.array())).array().sin();
               Array sinX = this->mSineA(0)*(this->mSineN(0)*(Math::PI/2)*(1+xGrid.array())).array().sin();
               MHDFloat yVal = std::cos(yGrid(iY));
               //MHDFloat yVal = 1.0;

               rNLComp.setProfile(sinZ*sinX(iX)*yVal,iY,iX);
            }
         }
      } else if(this->mTypeId == SINECOSINE)
      {
         int nX = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();
         int nY = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM2D,Dimensions::Space::PHYSICAL);
         int nZ = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM3D,Dimensions::Space::PHYSICAL);

         Array xGrid = Transform::TransformSelector<Dimensions::Transform::TRA1D>::Type::generateGrid(nX);
         Array yGrid = Transform::TransformSelector<Dimensions::Transform::TRA2D>::Type::generateGrid(nY);
         Array zGrid = Transform::TransformSelector<Dimensions::Transform::TRA3D>::Type::generateGrid(nZ);

         for(int iX = 0; iX < nX; ++iX)
         {
            for(int iY = 0; iY < nY; ++iY)
            {
               Array sinZ = this->mSineA(1)*(this->mSineN(1)*(Math::PI/2)*(1+zGrid.array())).array().cos();
               Array sinX = this->mSineA(0)*(this->mSineN(0)*(Math::PI/2)*(1+xGrid.array())).array().sin();
               MHDFloat yVal = std::cos(yGrid(iY));
               //MHDFloat yVal = 1.0;

               rNLComp.setProfile(sinZ*sinX(iX)*yVal,iY,iX);
            }
         }
      } else
      {
         throw Exception("Unknown exact state");
      }
   }

    Datatypes::SpectralScalarType::PointType ExactScalarState::sourceTerm(FieldComponents::Spectral::Id compId, const int i, const int j, const int k) const
    {
      // Assert on scalar component is used
      assert(compId == FieldComponents::Spectral::SCALAR);

      return Datatypes::SpectralScalarType::PointType(0,0);
    }

   void ExactScalarState::setRequirements()
   {
      // Add unknown to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(this->name(), FieldRequirement(true, true, true, false));
   }

}
}
