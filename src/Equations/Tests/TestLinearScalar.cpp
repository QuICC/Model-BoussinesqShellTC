/** 
 * @file TestLinearScalar.cpp
 * @brief Source of the implementation of the scalar linear test equation
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
#include "Equations/Tests/TestLinearScalar.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   TestLinearScalar::TestLinearScalar(SharedEquationParameters spEqParams)
      : IScalarEquation(spEqParams)
   {
   }

   TestLinearScalar::~TestLinearScalar()
   {
   }

   void TestLinearScalar::setIdentity(const PhysicalNames::Id name)
   {
      // Set the name
      this->setName(name);

      // Set the variable requirements
      this->setRequirements();
   }

   void TestLinearScalar::setCoupling()
   {
      this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::PROGNOSTIC, 0, false, false, false);
   }

   void TestLinearScalar::setRequirements()
   {
      // Add requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(this->name(), FieldRequirement(true, true, true, true));
   }

}
}
