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

   VectorFieldVisualizer::VectorFieldVisualizer(const std::string& pyName, SharedEquationParameters spEqParams)
      : IVectorEquation(pyName, spEqParams), mViewField(true), mViewGradient(false)
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
      for(SpectralComponent_iterator it = specRange.first; it != specRange.second; ++it)
      {
         // Initialise coupling information
         std::pair<std::map<FieldComponents::Spectral::Id, CouplingInformation>::iterator,bool> infoIt;
         infoIt = this->mCouplingInfos.insert(std::make_pair(*it,CouplingInformation()));
         SpectralFieldId eqId = std::make_pair(this->name(), *it);

         // General setup: first complex solver, complex solver, start from m = 0
         infoIt.first->second.setGeneral(CouplingInformation::TRIVIAL, false, 0);

         // Set nonlinear flags: NO nonlinear term, NO quasi-inverse
         infoIt.first->second.setNonlinear(false, false);

         // Set source flags: has source term
         infoIt.first->second.setSource(false);

         // Set index type: use MODE
         infoIt.first->second.setIndexType(CouplingInformation::MODE);

         // Equation is coupled to itself
         infoIt.first->second.addImplicitField(eqId.first, *it);

         // Set mininal matrix coupling
         int nMat = 0;
         ArrayI blockNs;
         ArrayI rhsCols;
         EigenSelector::makeMinimalCoupling(this->unknown().dom(0).spRes(), nMat, blockNs, rhsCols);
         infoIt.first->second.setSizes(nMat, blockNs, rhsCols); 

         // Sort implicit fields
         infoIt.first->second.sortImplicitFields(eqId.first, *it);
      }
   }

   void VectorFieldVisualizer::setRequirements()
   {
      // Add unknown to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(this->name(), FieldRequirement(false, true, this->mViewField, this->mViewGradient));
   }

}
}
