/** 
 * @file ScalarFieldVisualizer.cpp
 * @brief Source of the implementation of the basic scalar field visualizer
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
#include "Generator/Visualizers/ScalarFieldVisualizer.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "TypeSelectors/TransformSelector.hpp"
#include "TypeSelectors/EquationEigenSelector.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   ScalarFieldVisualizer::ScalarFieldVisualizer(SharedEquationParameters spEqParams)
      : IScalarEquation(spEqParams), mViewField(true), mViewGradient(false), mViewGradient2(false)
   {
   }

   ScalarFieldVisualizer::~ScalarFieldVisualizer()
   {
   }

   void ScalarFieldVisualizer::setIdentity(const PhysicalNames::Id name)
   {
      // Set the name
      this->setName(name);

      // Set the variable requirements
      this->setRequirements();
   }

   void ScalarFieldVisualizer::setFields(const bool viewField, const bool viewGradient, const bool viewGradient2)
   {
      this->mViewField = viewField;

      this->mViewGradient = viewGradient;

      this->mViewGradient2 = viewGradient2;
   }

   void ScalarFieldVisualizer::setCoupling()
   {
      this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::WRAPPER, 0, true, false);
   }

   void ScalarFieldVisualizer::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id compId) const
   {
   }

   void ScalarFieldVisualizer::setRequirements()
   {
      // Set solver timing
      this->setSolveTiming(SolveTiming::AFTER);

      // Add unknown to requirements: is scalar?, need spectral?, need physical?, need diff?, need curl?, need diff2?
      this->mRequirements.addField(this->name(), FieldRequirement(true, true, this->mViewField, this->mViewGradient, false, this->mViewGradient2));
   }

}
}
