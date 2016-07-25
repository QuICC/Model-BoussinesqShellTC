/** 
 * @file ScalarFieldTrivialVisualizer.cpp
 * @brief Source of the implementation of the basic scalar field visualizer with trivial solver
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
#include "Generator/Visualizers/ScalarFieldTrivialVisualizer.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "TypeSelectors/TransformSelector.hpp"
#include "Enums/FieldIdsTools.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   ScalarFieldTrivialVisualizer::ScalarFieldTrivialVisualizer(SharedEquationParameters spEqParams)
      : IScalarEquation(spEqParams), mViewField(true), mViewGradient(false), mViewGradient2(false)
   {
   }

   ScalarFieldTrivialVisualizer::~ScalarFieldTrivialVisualizer()
   {
   }

   void ScalarFieldTrivialVisualizer::setIdentity(const PhysicalNames::Id name)
   {
      // Set the name
      this->setName(name);

      // Set the variable requirements
      this->setRequirements();
   }

   void ScalarFieldTrivialVisualizer::setFields(const bool viewField, const bool viewGradient, const bool viewGradient2)
   {
      this->mViewField = viewField;

      this->mViewGradient = viewGradient;

      this->mViewGradient2 = viewGradient2;
   }

   void ScalarFieldTrivialVisualizer::setCoupling()
   {
      this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::TRIVIAL, 0, true, false);
   }

   void ScalarFieldTrivialVisualizer::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id compId) const
   {
   }

   void ScalarFieldTrivialVisualizer::setRequirements()
   {
      // Set solver timing
      this->setSolveTiming(SolveTiming::AFTER);

      // Add unknown to requirements: is scalar?, need spectral?, need physical?, need diff?, need curl?, need diff2?
      this->mRequirements.addField(this->name(), FieldRequirement(true, true, this->mViewField, this->mViewGradient, false, this->mViewGradient2));
   }

}
}
