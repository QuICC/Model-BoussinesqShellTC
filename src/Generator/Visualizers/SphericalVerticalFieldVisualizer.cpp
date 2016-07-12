/** 
 * @file SphericalVerticalFieldVisualizer.cpp
 * @brief Source of the implementation of the spherical vertical component field visualizer
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
#include "Generator/Visualizers/SphericalVerticalFieldVisualizer.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "TypeSelectors/TransformSelector.hpp"
#include "TypeSelectors/EquationEigenSelector.hpp"
#include "PhysicalOperators/SphericalZComponent.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   SphericalVerticalFieldVisualizer::SphericalVerticalFieldVisualizer(SharedEquationParameters spEqParams)
      : IScalarEquation(spEqParams), mFieldType(FieldType::VECTOR)
   {
   }

   SphericalVerticalFieldVisualizer::~SphericalVerticalFieldVisualizer()
   {
   }

   void SphericalVerticalFieldVisualizer::setIdentity(const PhysicalNames::Id vertName, const PhysicalNames::Id fieldName)
   {
      // Set the name
      this->setName(vertName);

      // Store field name
      this->mFieldName = fieldName;

      // Set the variable requirements
      this->setRequirements();
   }

   void SphericalVerticalFieldVisualizer::setFieldType(const FieldType::Id type)
   {
      if(type == FieldType::GRADIENT)
      {
         throw Exception("Z Component of gradient not implemented yet!");
      }

      this->mFieldType = type;
   }

   void SphericalVerticalFieldVisualizer::setCoupling()
   {
      this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::WRAPPER, 0, true, false);

      // Create cos(theta) and sin(theta) data for Coriolis term
      int nTh = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM2D,Dimensions::Space::PHYSICAL);
      Array thGrid = Transform::TransformSelector<Dimensions::Transform::TRA2D>::Type::generateGrid(nTh);
      this->mCosTheta = thGrid.array().cos();
      this->mSinTheta = thGrid.array().sin();
   }

   void SphericalVerticalFieldVisualizer::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id compId) const
   {
      if(this->mFieldType == FieldType::VECTOR)
      {
         Physical::SphericalZComponent::set(rNLComp, this->unknown().dom(0).spRes(), this->mCosTheta, this->mSinTheta, this->vector(mFieldName).dom(0).phys(), 1.0);
      } else if(this->mFieldType == FieldType::GRADIENT)
      {
         throw Exception("Z Component of gradient not implemented yet!");

         //Physical::SphericalZComponent::add(rNLComp, this->unknown().dom(0).spRes(), this->mCosTheta, this->mSinTheta, this->vector(mFieldName).dom(0).grad());
      } else if(this->mFieldType == FieldType::CURL)
      {
         Physical::SphericalZComponent::set(rNLComp, this->unknown().dom(0).spRes(), this->mCosTheta, this->mSinTheta, this->vector(mFieldName).dom(0).curl(), 1.0);
      }
   }

   void SphericalVerticalFieldVisualizer::useNonlinear(const Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id compId)
   {  
      this->rUnknown().rDom(0).rPhys().rData() = rNLComp.data();
   }

   void SphericalVerticalFieldVisualizer::setRequirements()
   {
      // Set solver timing
      this->setSolveTiming(SolveTiming::AFTER);

      // Add unknown to requirements: is scalar?, need spectral?, need physical?, need diff?(, needCurl)
      this->mRequirements.addField(this->name(), FieldRequirement(true, true, true, false, false, false));

      // Add base vector field to requirements: is scalar?, need spectral?, need physical?, need diff?(, needCurl)
      this->mRequirements.addField(this->mFieldName, FieldRequirement(false, true, (this->mFieldType == FieldType::VECTOR), (this->mFieldType == FieldType::GRADIENT), (this->mFieldType == FieldType::CURL)));
   }

}
}
