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
#include "TypeSelectors/SpectralSelector.hpp"
#include "TypeSelectors/TransformSelector.hpp"
#include "TypeSelectors/EquationToolsSelector.hpp"

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
      // Initialise coupling information
      std::pair<std::map<FieldComponents::Spectral::Id, CouplingInformation>::iterator,bool> infoIt;
      infoIt = this->mCouplingInfos.insert(std::make_pair(FieldComponents::Spectral::SCALAR,CouplingInformation()));
      SpectralFieldId eqId = std::make_pair(this->name(), FieldComponents::Spectral::SCALAR);

      // General setup: first complex solver, complex solver, start from m = 0
      infoIt.first->second.setGeneral(CouplingInformation::TRIVIAL, false, 0);

      // Set nonlinear flags: NO nonlinear term, NO quasi-inverse
      infoIt.first->second.setNonlinear(false, false);

      // Set source flags: has source term
      infoIt.first->second.setSource(false);

      // Equation is coupled to itself
      infoIt.first->second.addImplicitField(eqId.first, FieldComponents::Spectral::SCALAR);

      // Set mininal matrix coupling
      int nMat;
      ArrayI blockNs;
      ArrayI rhsCols;
      EquationToolsType::makeMinimalCoupling(this->unknown().dom(0).spRes(), nMat, blockNs, rhsCols);
      infoIt.first->second.setSizes(nMat, blockNs, rhsCols); 

      // Sort implicit fields
      infoIt.first->second.sortImplicitFields(eqId.first, FieldComponents::Spectral::SCALAR);
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
