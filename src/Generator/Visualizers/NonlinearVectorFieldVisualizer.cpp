/** 
 * @file NonlinearVectorFieldVisualizer.cpp
 * @brief Source of the implementation of a nonlinear vector field visualizer
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
#include "Generator/Visualizers/NonlinearVectorFieldVisualizer.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "TypeSelectors/TransformSelector.hpp"
#include "PhysicalOperators/Cross.hpp"

namespace QuICC {

namespace Equations {

   NonlinearVectorFieldVisualizer::NonlinearVectorFieldVisualizer(SharedEquationParameters spEqParams)
      : IVectorEquation(spEqParams)
   {
   }

   NonlinearVectorFieldVisualizer::~NonlinearVectorFieldVisualizer()
   {
   }

   void NonlinearVectorFieldVisualizer::setIdentity(const PhysicalNames::Id name)
   {
      // Set the name
      this->setName(name);

      // Set the variable requirements
      this->setRequirements();
   }

   void NonlinearVectorFieldVisualizer::setNonlinearType(const NonlinearVectorVisualizerIds::Id type)
   {
      this->mNonlinearType = type;
   }

   void NonlinearVectorFieldVisualizer::setCoupling()
   {
      if(FieldComponents::Spectral::ONE != FieldComponents::Spectral::NOTUSED)
      {
         this->defineCoupling(FieldComponents::Spectral::ONE, CouplingInformation::WRAPPER, 0, true, false);
      }

      if(FieldComponents::Spectral::TWO != FieldComponents::Spectral::NOTUSED)
      {
         this->defineCoupling(FieldComponents::Spectral::TWO, CouplingInformation::WRAPPER, 0, true, false);
      }

      if(FieldComponents::Spectral::THREE != FieldComponents::Spectral::NOTUSED)
      {
         this->defineCoupling(FieldComponents::Spectral::THREE, CouplingInformation::WRAPPER, 0, true, false);
      }
   }

   void NonlinearVectorFieldVisualizer::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id compId) const
   {
      if(this->mNonlinearType == NonlinearVectorVisualizerIds::CYLINDER_TORPOL_ADVECTION)
      {
         switch(compId)
         {
            case(FieldComponents::Physical::R):
               Physical::Cross<FieldComponents::Physical::THETA,FieldComponents::Physical::Z>::set(rNLComp, this->vector(this->name()).dom(0).curl(), this->vector(this->name()).dom(0).phys(), 1.0);
               break;
            case(FieldComponents::Physical::THETA):
               Physical::Cross<FieldComponents::Physical::Z,FieldComponents::Physical::R>::set(rNLComp, this->vector(this->name()).dom(0).curl(), this->vector(this->name()).dom(0).phys(), 1.0);
               break;
            case(FieldComponents::Physical::Z):
               Physical::Cross<FieldComponents::Physical::R,FieldComponents::Physical::THETA>::set(rNLComp, this->vector(this->name()).dom(0).curl(), this->vector(this->name()).dom(0).phys(), 1.0);
               break;
            default:
               assert(false);
               break;
         }
      }
   }

   void NonlinearVectorFieldVisualizer::useNonlinear(const Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id compId)
   {  
      this->rUnknown().rDom(0).rPhys().rComp(compId).rData() = rNLComp.data();
   }

   void NonlinearVectorFieldVisualizer::setRequirements()
   {
      // Set solver timing
      this->setSolveTiming(SolveTiming::AFTER);

      if(this->mNonlinearType == NonlinearVectorVisualizerIds::CYLINDER_TORPOL_ADVECTION)
      {
         // Add unknown to requirements: is scalar?, need spectral?, need physical?, need diff?(, needCurl)
         this->mRequirements.addField(this->name(), FieldRequirement(false, true, true, false, true));
      }
   }

   void NonlinearVectorFieldVisualizer::setNLComponents()
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
