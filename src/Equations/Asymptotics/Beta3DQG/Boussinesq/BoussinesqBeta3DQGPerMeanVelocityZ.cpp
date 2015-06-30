/** 
 * @file BoussinesqBeta3DQGPerMeanVelocityZ.cpp
 * @brief Source of the implementation of the mean Z velocity in the periodic Beta 3DQG model
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
#include "Equations/Asymptotics/Beta3DQG/Boussinesq/BoussinesqBeta3DQGPerMeanVelocityZ.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "Enums/NonDimensional.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   BoussinesqBeta3DQGPerMeanVelocityZ::BoussinesqBeta3DQGPerMeanVelocityZ(SharedEquationParameters spEqParams)
      : IScalarEquation(spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   BoussinesqBeta3DQGPerMeanVelocityZ::~BoussinesqBeta3DQGPerMeanVelocityZ()
   {
   }

   void BoussinesqBeta3DQGPerMeanVelocityZ::setCoupling()
   {
      this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::TRIVIAL, 0, false, false);
   }

   void BoussinesqBeta3DQGPerMeanVelocityZ::setRequirements()
   {
      // Set temperatur as equation unknown
      this->setName(PhysicalNames::MEAN_VELOCITYZ);

      // Set solver timing
      this->setSolveTiming(SolveTiming::AFTER);

      // Add mean Y velocity to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::MEAN_VELOCITYZ, FieldRequirement(true, true, true, false));

      // Add streamfunction to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::VELOCITYZ, FieldRequirement(true, true, false, false));
   }

}
}
