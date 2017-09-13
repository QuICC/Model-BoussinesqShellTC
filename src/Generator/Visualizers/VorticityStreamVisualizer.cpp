/** 
 * @file VorticityStreamVisualizer.cpp
 * @brief Source of the implementation of the streamfunction to vorticity visualizer
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Configuration includes
//

// System includes
//

// External includes

// Class include
//
#include "Generator/Visualizers/VorticityStreamVisualizer.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "Enums/NonDimensional.hpp"

namespace QuICC {

namespace Equations {

   VorticityStreamVisualizer::VorticityStreamVisualizer(SharedEquationParameters spEqParams)
      : IScalarEquation(spEqParams), mViewField(true), mViewGradient(false)
   {
   }

   VorticityStreamVisualizer::~VorticityStreamVisualizer()
   {
   }

   void VorticityStreamVisualizer::setFields(const bool viewField, const bool viewGradient)
   {
      this->mViewField = viewField;

      this->mViewGradient = viewGradient;

      // Set the variable requirements
      this->setRequirements();
   }

   void VorticityStreamVisualizer::setCoupling()
   {
      this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::WRAPPER, 0, false, false, false, false);
   }

   void VorticityStreamVisualizer::setRequirements()
   {
      // Set streamfunction as equation unknown
      this->setName(PhysicalNames::VORTICITY);

      // Set solver timing
      this->setSolveTiming(SolveTiming::AFTER);

      // Set vorticity requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::VORTICITY, FieldRequirement(true, true, this->mViewField, this->mViewGradient));

      // Add streamfunction requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::STREAMFUNCTION, FieldRequirement(true, true, false, false));
   }

}
}
