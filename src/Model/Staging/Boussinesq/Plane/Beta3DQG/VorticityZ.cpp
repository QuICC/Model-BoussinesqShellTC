/** 
 * @file BoussinesqBeta3DQGPerVorticityZ.cpp
 * @brief Source of the implementation of the vertical vorticity computation in the periodic Beta 3DQG model
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Configuration includes
//

// System includes
//

// External includes

// Class include
//
#include "Equations/Asymptotics/Beta3DQG/Boussinesq/BoussinesqBeta3DQGPerVorticityZ.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "Enums/NonDimensional.hpp"

namespace QuICC {

namespace Equations {

   BoussinesqBeta3DQGPerVorticityZ::BoussinesqBeta3DQGPerVorticityZ(SharedEquationParameters spEqParams)
      : IScalarEquation(spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   BoussinesqBeta3DQGPerVorticityZ::~BoussinesqBeta3DQGPerVorticityZ()
   {
   }

   void BoussinesqBeta3DQGPerVorticityZ::setCoupling()
   {
      this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::TRIVIAL, 1, false, false);
   }

   void BoussinesqBeta3DQGPerVorticityZ::setRequirements()
   {
      // Set streamfunction as equation unknown
      this->setName(PhysicalNames::VORTICITYZ);

      // Set solver timing
      this->setSolveTiming(SolveTiming::AFTER);

      // Add streamfunction requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::VORTICITYZ, FieldRequirement(true, true, false, false));
   }

}
}
