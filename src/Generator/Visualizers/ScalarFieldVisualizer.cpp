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
#include "TypeSelectors/SpectralOperatorSelector.hpp"
#include "TypeSelectors/TransformSelector.hpp"
#include "TypeSelectors/EquationEigenSelector.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   ScalarFieldVisualizer::ScalarFieldVisualizer(const std::string& pyName, SharedEquationParameters spEqParams)
      : IScalarEquation(pyName, spEqParams), mViewField(true), mViewGradient(false)
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

   void ScalarFieldVisualizer::setFields(const bool viewField, const bool viewGradient)
   {
      this->mViewField = viewField;

      this->mViewGradient = viewGradient;
   }

   void ScalarFieldVisualizer::setCoupling()
   {
      this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::WRAPPER, 0, false, false, false);
   }

   void ScalarFieldVisualizer::setRequirements()
   {
      // Add unknown to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(this->name(), FieldRequirement(true, true, this->mViewField, this->mViewGradient));
   }

}
}
