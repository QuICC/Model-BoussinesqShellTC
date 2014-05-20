/** 
 * @file BoussinesqBeta3DQGVorticity.cpp
 * @brief Source of the implementation of the vorticity computation in the Beta 3DQG model
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Configuration includes
//

// System includes
//

// External includes

// Class include
//
#include "Equations/Asymptotics/Beta3DQG/Boussinesq/BoussinesqBeta3DQGVorticity.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "Enums/NonDimensional.hpp"
#include "SpectralOperators/Tools/SpectralBoxTools.hpp"
#include "TypeSelectors/SpectralOperatorSelector.hpp"
#include "TypeSelectors/EquationEigenSelector.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   BoussinesqBeta3DQGVorticity::BoussinesqBeta3DQGVorticity(const std::string& pyName, SharedEquationParameters spEqParams)
      : IScalarEquation(pyName,spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   BoussinesqBeta3DQGVorticity::~BoussinesqBeta3DQGVorticity()
   {
   }

   void BoussinesqBeta3DQGVorticity::setCoupling()
   {
      this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::TRIVIAL, 0, false, false, false);
   }

   void BoussinesqBeta3DQGVorticity::setRequirements()
   {
      // Set streamfunction as equation unknown
      this->setName(PhysicalNames::VORTICITY);

      // Set vorticity requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::VORTICITY, FieldRequirement(true, true, true, true));

      // Add streamfunction requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::STREAMFUNCTION, FieldRequirement(true, true, false, false));
   }

}
}
