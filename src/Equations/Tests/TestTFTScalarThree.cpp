/** \file TestTFTScalarThree.cpp
 *  \brief Source of the implementation of the third scalar test equation for the TFT scheme
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "Equations/Tests/TestTFTScalarThree.hpp"

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

   TestTFTScalarThree::TestTFTScalarThree(SharedEquationParameters spEqParams)
      : ITestTFTScalarEquation(spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   TestTFTScalarThree::~TestTFTScalarThree()
   {
   }

   void TestTFTScalarThree::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const
   {
      // Assert on scalar component is used
      assert(id == FieldComponents::Physical::SCALAR);

      rNLComp.rData().setConstant(0.0);
   }

   void TestTFTScalarThree::setRequirements()
   {
      // Set temperature as equation unknown
      this->setName(PhysicalNames::VELOCITYZ);

      // Add temperature to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::VELOCITYZ, FieldRequirement(true, true, true, true));
   }

}
}
