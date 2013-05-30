/** \file TestTFTScalarTwo.cpp
 *  \brief Source of the implementation of the second scalar test equation for the TFT scheme
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "Equations/Tests/TestTFTScalarTwo.hpp"

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

   TestTFTScalarTwo::TestTFTScalarTwo(SharedEquationParameters spEqParams)
      : ITestTFTScalarEquation(spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   TestTFTScalarTwo::~TestTFTScalarTwo()
   {
   }

   void TestTFTScalarTwo::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp) const
   {
      rNLComp.rData().setConstant(0.0);
   }

   void TestTFTScalarTwo::setRequirements()
   {
      // Set temperature as equation unknown
      this->setName(PhysicalNames::STREAMFUNCTION);

      // Add temperature to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::STREAMFUNCTION, FieldRequirement(true, true, true, true));
   }

}
}
