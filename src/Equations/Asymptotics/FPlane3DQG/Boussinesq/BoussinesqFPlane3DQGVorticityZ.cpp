/** 
 * @file BoussinesqFPlane3DQGVorticityZ.cpp
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
#include "Equations/Asymptotics/FPlane3DQG/Boussinesq/BoussinesqFPlane3DQGVorticityZ.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "Enums/NonDimensional.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   BoussinesqFPlane3DQGVorticityZ::BoussinesqFPlane3DQGVorticityZ(SharedEquationParameters spEqParams)
      : IScalarEquation(spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   BoussinesqFPlane3DQGVorticityZ::~BoussinesqFPlane3DQGVorticityZ()
   {
   }

   void BoussinesqFPlane3DQGVorticityZ::setCoupling()
   {
      this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::TRIVIAL, 1, false, false, false);
      this->setExplicitTiming(FieldComponents::Spectral::SCALAR, ExplicitTiming::LINEAR);
   }

   void BoussinesqFPlane3DQGVorticityZ::setRequirements()
   {
      // Set streamfunction as equation unknown
      this->setName(PhysicalNames::VORTICITYZ);

      // Set solver timing
      this->setSolveTiming(SolveTiming::AFTER);

      // Set non orthogonal vertical vorticity requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::VORTICITYZ, FieldRequirement(true, true, true, true));
   }

}
}
