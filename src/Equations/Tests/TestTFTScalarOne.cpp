/** \file TestTFTScalarOne.cpp
 *  \brief Source of the implementation of the first scalar test equation for the TFT scheme
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "Equations/Tests/TestTFTScalarOne.hpp"

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

   TestTFTScalarOne::TestTFTScalarOne(SharedEquationParameters spEqParams)
      : ITestTFTScalarEquation(spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   TestTFTScalarOne::~TestTFTScalarOne()
   {
   }

   void TestTFTScalarOne::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const
   {
      // Assert on scalar component is used
      assert(id == FieldComponents::Physical::SCALAR);

      rNLComp.rData().setConstant(0.0);
   }

   void TestTFTScalarOne::setRequirements()
   {
      // Set temperature as equation unknown
      this->setName(PhysicalNames::TEMPERATURE);

      // Add temperature to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::TEMPERATURE, FieldRequirement(true, true, true, true));
   }

}
}
