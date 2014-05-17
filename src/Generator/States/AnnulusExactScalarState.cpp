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
#include "TypeSelectors/SpectralOperatorSelector.hpp"
#include "TypeSelectors/TransformSelector.hpp"
#include "TypeSelectors/EquationEigenSelector.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   AnnulusExactScalarState::AnnulusExactScalarState(const std::string& pyName, SharedEquationParameters spEqParams)
      : IScalarEquation(pyName,spEqParams), mTypeId(CONSTANT)
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

   void AnnulusExactScalarState::setStateType(const AnnulusExactScalarState::StateTypeId id)
   {
      this->mTypeId = id;
   }

   void AnnulusExactScalarState::setCoupling()
   {
      this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::TRIVIAL, 0, true, false, false);
   }

   void AnnulusExactScalarState::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id compId) const
   {
      // Assert on scalar component is used
      assert(compId == FieldComponents::Physical::SCALAR);

      if(this->mTypeId == CONSTANT)
      {
         int nR = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();
         int nTh = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM2D,Dimensions::Space::PHYSICAL);
         int nZ = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM3D,Dimensions::Space::PHYSICAL);

         Array rGrid = Transform::TransformSelector<Dimensions::Transform::TRA1D>::Type::generateGrid(nR, 1, 0.35);
         Array thGrid = Transform::TransformSelector<Dimensions::Transform::TRA2D>::Type::generateGrid(nTh);
         Array zGrid = Transform::TransformSelector<Dimensions::Transform::TRA3D>::Type::generateGrid(nZ);

         for(int iR = 0; iR < nR; ++iR)
         {
            for(int iTh = 0; iTh < nTh; ++iTh)
            {
               Array funcZ = 2.*(5.*(Math::PI/2)*(1+zGrid.array())).array().sin();
               Array funcR = 1 + rGrid.array() + rGrid.array().pow(3);
               MHDFloat funcTh = std::cos(thGrid(iTh));

               rNLComp.setProfile(funcZ*funcR(iR)*funcTh,iTh,iR);
            }
         }
      } else
      {
         throw Exception("Unknown exact state");
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
      // Add unknown to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(this->name(), FieldRequirement(true, true, false, false));
   }

}
}
