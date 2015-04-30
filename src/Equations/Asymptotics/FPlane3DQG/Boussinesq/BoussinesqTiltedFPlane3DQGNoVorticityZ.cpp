/** 
 * @file BoussinesqTiltedFPlane3DQGNoVorticityZ.cpp
 * @brief Source of the implementation of the non orthogonal vertical vorticity computation in the tilted F-plane 3DQG model
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Configuration includes
//

// System includes
//

// External includes

// Class include
//
#include "Equations/Asymptotics/FPlane3DQG/Boussinesq/BoussinesqTiltedFPlane3DQGNoVorticityZ.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "Enums/NonDimensional.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   BoussinesqTiltedFPlane3DQGNoVorticityZ::BoussinesqTiltedFPlane3DQGNoVorticityZ(SharedEquationParameters spEqParams)
      : IScalarEquation(spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   BoussinesqTiltedFPlane3DQGNoVorticityZ::~BoussinesqTiltedFPlane3DQGNoVorticityZ()
   {
   }

   void BoussinesqTiltedFPlane3DQGNoVorticityZ::setCoupling()
   {
      this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::TRIVIAL, 1, false, false);
   }

   void BoussinesqTiltedFPlane3DQGNoVorticityZ::setRequirements()
   {
      // Set streamfunction as equation unknown
      this->setName(PhysicalNames::NO_VORTICITYZ);

      // Set solver timing
      this->setSolveTiming(SolveTiming::AFTER);

      // Set non orthogonal vertical vorticity requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::NO_VORTICITYZ, FieldRequirement(true, true, true, true));
   }

}
}
