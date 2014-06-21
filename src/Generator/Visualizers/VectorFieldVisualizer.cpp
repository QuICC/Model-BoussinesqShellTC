/** 
 * @file VectorFieldVisualizer.cpp
 * @brief Source of the implementation of the basic vector field visualizer
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
#include "Generator/Visualizers/VectorFieldVisualizer.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "TypeSelectors/SpectralOperatorSelector.hpp"
#include "TypeSelectors/TransformSelector.hpp"
#include "TypeSelectors/EquationEigenSelector.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   VectorFieldVisualizer::VectorFieldVisualizer(SharedEquationParameters spEqParams)
      : IVectorEquation(spEqParams), mViewField(true), mViewGradient(false)
   {
   }

   VectorFieldVisualizer::~VectorFieldVisualizer()
   {
   }

   void VectorFieldVisualizer::setIdentity(const PhysicalNames::Id name)
   {
      // Set the name
      this->setName(name);

      // Setup toroidal/poloidal components
      this->mSpectralIds.push_back(FieldComponents::Spectral::ONE);
      this->mSpectralIds.push_back(FieldComponents::Spectral::TWO);

      // Set the variable requirements
      this->setRequirements();
   }

   void VectorFieldVisualizer::setFields(const bool viewField, const bool viewGradient)
   {
      this->mViewField = viewField;

      this->mViewGradient = viewGradient;
   }

   void VectorFieldVisualizer::setCoupling()
   {
      SpectralComponent_range specRange = this->spectralRange();
      SpectralComponent_iterator specIt;
      for(specIt = specRange.first; specIt != specRange.second; ++specIt)
      {
         this->defineCoupling(*specIt, CouplingInformation::WRAPPER, 0, false, false, false);
      }
   }

   void VectorFieldVisualizer::setRequirements()
   {
      // Add unknown to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(this->name(), FieldRequirement(false, true, this->mViewField, this->mViewGradient));
   }

}
}
