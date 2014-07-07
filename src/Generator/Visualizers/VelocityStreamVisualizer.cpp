/** 
 * @file VelocityStreamVisualizer.cpp
 * @brief Source of the implementation of the streamfunction + axial velocity to velocity field visualizer
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
#include "Generator/Visualizers/VelocityStreamVisualizer.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "TypeSelectors/TransformSelector.hpp"
#include "TypeSelectors/EquationEigenSelector.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   VelocityStreamVisualizer::VelocityStreamVisualizer(SharedEquationParameters spEqParams)
      : IVectorEquation(spEqParams), mViewField(true), mViewGradient(false)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   VelocityStreamVisualizer::~VelocityStreamVisualizer()
   {
   }

   void VelocityStreamVisualizer::setIdentity(const PhysicalNames::Id name)
   {
   }

   void VelocityStreamVisualizer::setFields(const bool viewField, const bool viewGradient)
   {
      this->mViewField = viewField;

      this->mViewGradient = viewGradient;
   }

   void VelocityStreamVisualizer::setCoupling()
   {
      SpectralComponent_range specRange = this->spectralRange();
      SpectralComponent_iterator specIt;
      for(specIt = specRange.first; specIt != specRange.second; ++specIt)
      {
         this->defineCoupling(*specIt, CouplingInformation::WRAPPER, 0, false, false, false);
      }
   }

   void VelocityStreamVisualizer::setRequirements()
   {
      // Set the name
      this->setName(PhysicalNames::VELOCITY);

      // Add unknown to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(this->name(), FieldRequirement(false, false, this->mViewField, this->mViewGradient));
   }

}
}
