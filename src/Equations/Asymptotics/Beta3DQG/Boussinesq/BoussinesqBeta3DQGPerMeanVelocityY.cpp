/** 
 * @file BoussinesqBeta3DQGPerMeanVelocityY.cpp
 * @brief Source of the implementation of the mean Y velocity in the periodic Beta 3DQG model
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
#include "Equations/Asymptotics/Beta3DQG/Boussinesq/BoussinesqBeta3DQGPerMeanVelocityY.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "Enums/NonDimensional.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   BoussinesqBeta3DQGPerMeanVelocityY::BoussinesqBeta3DQGPerMeanVelocityY(SharedEquationParameters spEqParams)
      : IScalarEquation(spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   BoussinesqBeta3DQGPerMeanVelocityY::~BoussinesqBeta3DQGPerMeanVelocityY()
   {
   }

   void BoussinesqBeta3DQGPerMeanVelocityY::setCoupling()
   {
      this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::TRIVIAL, 0, false, false);
   }

   void BoussinesqBeta3DQGPerMeanVelocityY::setRequirements()
   {
      // Set temperatur as equation unknown
      this->setName(PhysicalNames::MEAN_VELOCITYY);

      // Set solver timing
      this->setSolveTiming(SolveTiming::AFTER);

      // Add mean Y velocity to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::MEAN_VELOCITYY, FieldRequirement(true, true, true, false));

      // Add streamfunction to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::STREAMFUNCTION, FieldRequirement(true, true, false, false));
   }

}
}
