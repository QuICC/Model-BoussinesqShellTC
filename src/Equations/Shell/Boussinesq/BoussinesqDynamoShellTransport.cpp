/** 
 * @file BoussinesqDynamoShellTransport.cpp
 * @brief Source of the implementation of the transport equation in the Boussinesq thermal convection dynamo in a spherical shell
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
#include "Equations/Shell/Boussinesq/BoussinesqDynamoShellTransport.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   BoussinesqDynamoShellTransport::BoussinesqDynamoShellTransport(SharedEquationParameters spEqParams)
      : IScalarEquation(spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   BoussinesqDynamoShellTransport::~BoussinesqDynamoShellTransport()
   {
   }

   void BoussinesqDynamoShellTransport::setCoupling()
   {
      this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::PROGNOSTIC, 0, false, false, false);
   }

   void BoussinesqDynamoShellTransport::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const
   {
      // Assert on scalar component is used
      assert(id == FieldComponents::Physical::SCALAR);

      throw Exception("Nonlinear term in spherical shell transport equation not done yet");
   }

   void BoussinesqDynamoShellTransport::setRequirements()
   {
      // Set temperatur as equation unknown
      this->setName(PhysicalNames::TEMPERATURE);

      // Set solver timing
      this->setSolveTiming(SolveTiming::PROGNOSTIC);

      // Add temperature to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::TEMPERATURE, FieldRequirement(true, true, false, false));
   }

}
}
