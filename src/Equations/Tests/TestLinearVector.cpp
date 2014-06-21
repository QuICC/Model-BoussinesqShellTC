/** 
 * @file TestLinearVector.cpp
 * @brief Source of the implementation of the vector linear test equation
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
#include "Equations/Tests/TestLinearVector.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   TestLinearVector::TestLinearVector(SharedEquationParameters spEqParams)
      : IVectorEquation(spEqParams)
   {
   }

   TestLinearVector::~TestLinearVector()
   {
   }

   void TestLinearVector::setIdentity(const PhysicalNames::Id name)
   {
      // Set the name
      this->setName(name);

      // Set the variable requirements
      this->setRequirements();
   }

   void TestLinearVector::setCoupling()
   {
      SpectralComponent_range specRange = this->spectralRange();
      SpectralComponent_iterator specIt;
      for(specIt = specRange.first; specIt != specRange.second; ++specIt)
      {
         this->defineCoupling(*specIt, CouplingInformation::PROGNOSTIC, 1, false, false, false);
      }
   }

   void TestLinearVector::setRequirements()
   {
      // Add requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(this->name(), FieldRequirement(false, true, true, true));
   }

}
}
