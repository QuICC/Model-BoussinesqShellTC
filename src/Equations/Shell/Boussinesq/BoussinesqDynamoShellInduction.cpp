/** 
 * @file BoussinesqDynamoShellInduction.cpp
 * @brief Source of the implementation of the vector induction equation in the Boussinesq thermal convection dynamo in a spherical shell model
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
#include "Equations/Shell/Boussinesq/BoussinesqDynamoShellInduction.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "Enums/NonDimensional.hpp"
#include "TypeSelectors/EquationEigenSelector.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   BoussinesqDynamoShellInduction::BoussinesqDynamoShellInduction(SharedEquationParameters spEqParams)
      : IVectorEquation(spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   BoussinesqDynamoShellInduction::~BoussinesqDynamoShellInduction()
   {
   }

   void BoussinesqDynamoShellInduction::setCoupling()
   {
      this->defineCoupling(FieldComponents::Spectral::TOR, CouplingInformation::PROGNOSTIC, 1, false, false, false);

      this->defineCoupling(FieldComponents::Spectral::POL, CouplingInformation::PROGNOSTIC, 1, false, false, false);
   }

   void BoussinesqDynamoShellInduction::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const
   {
      throw Exception("Nonlinear term in spherical shell toroidal/poloidal induction equation not done yet");
   }

   void BoussinesqDynamoShellInduction::setRequirements()
   {
      // Set velocity as equation unknown
      this->setName(PhysicalNames::MAGNETIC);

      // Set solver timing
      this->setSolveTiming(SolveTiming::PROGNOSTIC);

      // Add velocity to requirements: is scalar?, need spectral?, need physical?, need diff?(, need curl?)
      this->mRequirements.addField(PhysicalNames::MAGNETIC, FieldRequirement(false, true, false, false));
   }

}
}
