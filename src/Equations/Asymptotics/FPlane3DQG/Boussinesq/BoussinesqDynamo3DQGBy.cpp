/** 
 * @file BoussinesqDynamo3DQGBy.cpp
 * @brief Source of the implementation of the vertical vorticity computation in the F-plane 3DQG model
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Configuration includes
//

// System includes
//

// External includes

// Class include
//
#include "Equations/Asymptotics/FPlane3DQG/Boussinesq/BoussinesqDynamo3DQGBy.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "Enums/NonDimensional.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   BoussinesqDynamo3DQGBy::BoussinesqDynamo3DQGBy(SharedEquationParameters spEqParams)
      : IScalarEquation(spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   BoussinesqDynamo3DQGBy::~BoussinesqDynamo3DQGBy()
   {
   }

   void BoussinesqDynamo3DQGBy::setCoupling()
   {
      this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::PROGNOSTIC, 0, false, false);
   }

   void BoussinesqDynamo3DQGBy::setRequirements()
   {
      // Set streamfunction as equation unknown
      this->setName(PhysicalNames::BY);

      // Set solver timing
      this->setSolveTiming(SolveTiming::PROGNOSTIC);

      // Set non orthogonal vertical vorticity requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::BY, FieldRequirement(true, true, false, false));
      // Set non orthogonal vertical vorticity requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::EMFX, FieldRequirement(true, true, false, false));
   }

}
}
