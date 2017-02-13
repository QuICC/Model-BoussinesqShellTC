/** 
 * @file BoussinesqFPlaneNHBGEVorticityZ.cpp
 * @brief Source of the implementation of the vertical vorticity computation in the F-plane NHBGE model
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Configuration includes
//

// System includes
//

// External includes

// Class include
//
#include "Equations/Asymptotics/FPlaneNHBGE/Boussinesq/BoussinesqFPlaneNHBGEVorticityZ.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "Enums/NonDimensional.hpp"

namespace QuICC {

namespace Equations {

   BoussinesqFPlaneNHBGEVorticityZ::BoussinesqFPlaneNHBGEVorticityZ(SharedEquationParameters spEqParams)
      : IScalarEquation(spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   BoussinesqFPlaneNHBGEVorticityZ::~BoussinesqFPlaneNHBGEVorticityZ()
   {
   }

   void BoussinesqFPlaneNHBGEVorticityZ::setCoupling()
   {
      this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::TRIVIAL, 1, false, false);
   }

   void BoussinesqFPlaneNHBGEVorticityZ::setRequirements()
   {
      // Set streamfunction as equation unknown
      this->setName(PhysicalNames::VORTICITYZ);

      // Set solver timing
      this->setSolveTiming(SolveTiming::AFTER);

      // Set non orthogonal vertical vorticity requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::VORTICITYZ, FieldRequirement(true, true, false, false));
   }

}
}
