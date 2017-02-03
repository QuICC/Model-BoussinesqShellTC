/** 
 * @file BoussinesqDynamoCouetteShellInduction.cpp
 * @brief Source of the implementation of the vector induction equation in the Boussinesq spherical Couette dynamo in a spherical shell model
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
#include "Equations/Shell/Boussinesq/BoussinesqDynamoCouetteShellInduction.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "Enums/NonDimensional.hpp"
#include "PhysicalOperators/Cross.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   BoussinesqDynamoCouetteShellInduction::BoussinesqDynamoCouetteShellInduction(SharedEquationParameters spEqParams)
      : IVectorEquation(spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   BoussinesqDynamoCouetteShellInduction::~BoussinesqDynamoCouetteShellInduction()
   {
   }

   void BoussinesqDynamoCouetteShellInduction::setCoupling()
   {
      #ifdef GEOMHDISCC_SPATIALSCHEME_SLFL
         int start = 1;
      #else //if GEOMHDISCC_SPATIALSCHEME_SLFM
         int start = 0;
      #endif //GEOMHDISCC_SPATIALSCHEME_SLFL

      this->defineCoupling(FieldComponents::Spectral::TOR, CouplingInformation::PROGNOSTIC, start, true, false);

      this->defineCoupling(FieldComponents::Spectral::POL, CouplingInformation::PROGNOSTIC, start, true, false);
   }

   void BoussinesqDynamoCouetteShellInduction::setNLComponents()
   {
      this->addNLComponent(FieldComponents::Spectral::POL,0);

      this->addNLComponent(FieldComponents::Spectral::TOR,1);
   }

   void BoussinesqDynamoCouetteShellInduction::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id compId) const
   {
	  // Import the Rossby number for scaling the induction term
	   MHDFloat Ro = std::abs(this->eqParams().nd(NonDimensional::ROSSBY));
      ///
      /// Compute \f$\left(\vec u\wedge\vec B + B0\right)\f$
      ///
      switch(compId)
      {
         case(FieldComponents::Physical::R):
            Physical::Cross<FieldComponents::Physical::THETA,FieldComponents::Physical::PHI>::set(rNLComp, this->unknown().dom(0).phys(), this->vector(PhysicalNames::VELOCITY).dom(0).phys(), Ro);
            Physical::Cross<FieldComponents::Physical::THETA,FieldComponents::Physical::PHI>::add(rNLComp, this->vector(PhysicalNames::IMPOSED_MAGNETIC).dom(0).phys(), this->vector(PhysicalNames::VELOCITY).dom(0).phys(), Ro);
            break;
         case(FieldComponents::Physical::THETA):
            Physical::Cross<FieldComponents::Physical::PHI,FieldComponents::Physical::R>::set(rNLComp, this->unknown().dom(0).phys(), this->vector(PhysicalNames::VELOCITY).dom(0).phys(), Ro);
            Physical::Cross<FieldComponents::Physical::PHI,FieldComponents::Physical::R>::add(rNLComp, this->vector(PhysicalNames::IMPOSED_MAGNETIC).dom(0).phys(), this->vector(PhysicalNames::VELOCITY).dom(0).phys(), Ro);
            break;
         case(FieldComponents::Physical::PHI):
            Physical::Cross<FieldComponents::Physical::R,FieldComponents::Physical::THETA>::set(rNLComp, this->unknown().dom(0).phys(), this->vector(PhysicalNames::VELOCITY).dom(0).phys(), Ro);
            Physical::Cross<FieldComponents::Physical::R,FieldComponents::Physical::THETA>::add(rNLComp, this->vector(PhysicalNames::IMPOSED_MAGNETIC).dom(0).phys(), this->vector(PhysicalNames::VELOCITY).dom(0).phys(), Ro);
            break;
         default:
            assert(false);
            break;
      }
   }

   void BoussinesqDynamoCouetteShellInduction::setRequirements()
   {
      // Set velocity as equation unknown
      this->setName(PhysicalNames::MAGNETIC);

      // Set solver timing
      this->setSolveTiming(SolveTiming::PROGNOSTIC);

      // Add magnetic field to requirements: is scalar?, need spectral?, need physical?, need diff?(, need curl?)
      this->mRequirements.addField(PhysicalNames::MAGNETIC, FieldRequirement(false, true, true, false, false));

      // Add velocity to requirements: is scalar?, need spectral?, need physical?, need diff?(, need curl?)
      this->mRequirements.addField(PhysicalNames::VELOCITY, FieldRequirement(false, true, true, false, false));

      // Add imposed magnetic field to requirements: is scalar?, need spectral?, need physical?, need diff?(, need curl?)
      this->mImposedRequirements.addField(PhysicalNames::IMPOSED_MAGNETIC, FieldRequirement(false, true, true, false, false));
   }

}
}
