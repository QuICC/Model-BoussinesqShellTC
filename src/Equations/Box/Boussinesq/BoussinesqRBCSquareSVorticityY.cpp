/** 
 * @file BoussinesqRBCSquareSVorticityY.cpp
 * @brief Source of the implementation of the vorticity equation in Rayleigh-Benard convection in a square cavity (2D) (streamfunction formulation)
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
#include "Equations/Box/Boussinesq/BoussinesqRBCSquareSVorticityY.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "Enums/NonDimensional.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   BoussinesqRBCSquareSVorticityY::BoussinesqRBCSquareSVorticityY(SharedEquationParameters spEqParams)
      : IScalarEquation(spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   BoussinesqRBCSquareSVorticityY::~BoussinesqRBCSquareSVorticityY()
   {
   }

   void BoussinesqRBCSquareSVorticityY::setCoupling()
   {
      this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::DIAGNOSTIC, 0, false, false);
   }

   void BoussinesqRBCSquareSVorticityY::setRequirements()
   {
      // Set temperatur as equation unknown
      this->setName(PhysicalNames::VORTICITYY);

      // Set solver timing
      this->setSolveTiming(SolveTiming::AFTER);

      // Add temperature to requirements: is scalar?, need spectral?, need physical?, need diff?(, need curl?)
      this->mRequirements.addField(PhysicalNames::VORTICITYY, FieldRequirement(true, true, false, true));
   }

}
}
