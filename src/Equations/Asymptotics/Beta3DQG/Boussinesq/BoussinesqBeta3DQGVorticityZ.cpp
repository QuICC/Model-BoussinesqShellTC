/** 
 * @file BoussinesqBeta3DQGVorticityZ.cpp
 * @brief Source of the implementation of the vertical vorticity computation in the Beta 3DQG model
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Configuration includes
//

// System includes
//

// External includes

// Class include
//
#include "Equations/Asymptotics/Beta3DQG/Boussinesq/BoussinesqBeta3DQGVorticityZ.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "Enums/NonDimensional.hpp"
#include "TypeSelectors/EquationEigenSelector.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   BoussinesqBeta3DQGVorticityZ::BoussinesqBeta3DQGVorticityZ(SharedEquationParameters spEqParams)
      : IScalarEquation(spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   BoussinesqBeta3DQGVorticityZ::~BoussinesqBeta3DQGVorticityZ()
   {
   }

   void BoussinesqBeta3DQGVorticityZ::setCoupling()
   {
      this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::TRIVIAL, 0, false, false, false);
   }

   void BoussinesqBeta3DQGVorticityZ::setRequirements()
   {
      // Set streamfunction as equation unknown
      this->setName(PhysicalNames::VORTICITYZ);

      // Set vorticity requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::VORTICITYZ, FieldRequirement(true, true, true, true));

      // Add streamfunction requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::STREAMFUNCTION, FieldRequirement(true, true, false, false));
   }

}
}
