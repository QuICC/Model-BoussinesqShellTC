/** \file TestTFTTransport.cpp
 *  \brief Source of the implementation of the test equation for the TFT scheme
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "Equations/Tests/TestTFTTransport.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "PhysicalOperators/StreamAdvection.hpp"
#include "SpectralOperators/PeriodicOperator.hpp"
#include "TypeSelectors/SpectralSelector.hpp"
#include "Equations/Tests/TestTFTSystem.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   TestTFTTransport::TestTFTTransport(SharedEquationParameters spEqParams)
      : ITestTFTScalarEquation(spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   TestTFTTransport::~TestTFTTransport()
   {
   }

   void TestTFTTransport::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp) const
   {
   }

   void TestTFTTransport::setRequirements()
   {
      // Set temperature as equation unknown
      this->setName(PhysicalNames::TEMPERATURE);

      // Add temperature to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::TEMPERATURE, FieldRequirement(true, true, false, false));
   }

}
}
