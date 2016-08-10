/** 
 * @file VectorFieldTrivialVisualizer.cpp
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
#include "Generator/Visualizers/VectorFieldTrivialVisualizer.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "TypeSelectors/TransformSelector.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   VectorFieldTrivialVisualizer::VectorFieldTrivialVisualizer(SharedEquationParameters spEqParams)
      : IVectorEquation(spEqParams), mViewField(true), mViewGradient(false), mViewCurl(false)
   {
   }

   VectorFieldTrivialVisualizer::~VectorFieldTrivialVisualizer()
   {
   }

   void VectorFieldTrivialVisualizer::setIdentity(const PhysicalNames::Id name)
   {
      // Set the name
      this->setName(name);

      // Set the variable requirements
      this->setRequirements();
   }

   void VectorFieldTrivialVisualizer::setFields(const bool viewField, const bool viewGradient, const bool viewCurl)
   {
      this->mViewField = viewField;

      this->mViewGradient = viewGradient;

      this->mViewCurl = viewCurl;
   }

   void VectorFieldTrivialVisualizer::setCoupling()
   {
      if(FieldComponents::Spectral::ONE != FieldComponents::Spectral::NOTUSED)
      {
         this->defineCoupling(FieldComponents::Spectral::ONE, CouplingInformation::TRIVIAL, 0, true, false);
      }

      if(FieldComponents::Spectral::TWO != FieldComponents::Spectral::NOTUSED)
      {
         this->defineCoupling(FieldComponents::Spectral::TWO, CouplingInformation::TRIVIAL, 0, true, false);
      }

      if(FieldComponents::Spectral::THREE != FieldComponents::Spectral::NOTUSED)
      {
         this->defineCoupling(FieldComponents::Spectral::THREE, CouplingInformation::TRIVIAL, 0, true, false);
      }
   }

   void VectorFieldTrivialVisualizer::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id compId) const
   {
   }

   void VectorFieldTrivialVisualizer::useNonlinear(const Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id compId)
   {  
   }

   void VectorFieldTrivialVisualizer::setRequirements()
   {
      // Set solver timing
      this->setSolveTiming(SolveTiming::AFTER);

      // Add unknown to requirements: is scalar?, need spectral?, need physical?, need diff?(, needCurl)
      this->mRequirements.addField(this->name(), FieldRequirement(false, true, this->mViewField, this->mViewGradient, this->mViewCurl));
   }

   void VectorFieldTrivialVisualizer::setNLComponents()
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
