/** 
 * @file NonlinearScalarFieldVisualizer.cpp
 * @brief Source of the implementation of a nonlinearscalar field visualizer
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
#include "Generator/Visualizers/NonlinearScalarFieldVisualizer.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "TypeSelectors/TransformSelector.hpp"
#include "Enums/FieldIdsTools.hpp"
#include "PhysicalOperators/VelocityHeatAdvection.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   NonlinearScalarFieldVisualizer::NonlinearScalarFieldVisualizer(SharedEquationParameters spEqParams)
      : IScalarEquation(spEqParams)
   {
   }

   NonlinearScalarFieldVisualizer::~NonlinearScalarFieldVisualizer()
   {
   }

   void NonlinearScalarFieldVisualizer::setIdentity(const PhysicalNames::Id name)
   {
      // Set the name
      this->setName(name);

      // Set the variable requirements
      this->setRequirements();
   }

   void NonlinearScalarFieldVisualizer::setNonlinearType(const NonlinearScalarVisualizerIds::Id type)
   {
      this->mNonlinearType = type;
   }

   void NonlinearScalarFieldVisualizer::setCoupling()
   {
      this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::WRAPPER, 0, true, false);
   }

   void NonlinearScalarFieldVisualizer::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id compId) const
   {
      if(this->mNonlinearType == NonlinearScalarVisualizerIds::CYLINDER_HEAT_ADVECTION)
      {
         // Assert on scalar component is used
         assert(id == FieldComponents::Physical::SCALAR);                                                                                                        
         Physical::VelocityHeatAdvection<FieldComponents::Physical::R,FieldComponents::Physical::THETA,FieldComponents::Physical::Z>::set(rNLComp, this->vector(PhysicalNames::VELOCITY).dom(0).phys(), this->scalar(PhysicalNames::TEMPERATURE).dom(0).grad(), 1.0);
      }
   }

   void NonlinearScalarFieldVisualizer::useNonlinear(const Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id compId)
   {  
      this->rUnknown().rDom(0).rPhys().rData() = rNLComp.data();
   }

   void NonlinearScalarFieldVisualizer::setRequirements()
   {
      // Set solver timing
      this->setSolveTiming(SolveTiming::AFTER);

      // Add unknown to requirements: is scalar?, need spectral?, need physical?, need diff?, need curl?, need diff2?
      this->mRequirements.addField(this->name(), FieldRequirement(true, true, this->mViewField, this->mViewGradient, false, this->mViewGradient2));
   }

}
}
