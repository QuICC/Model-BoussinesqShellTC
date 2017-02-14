/** 
 * @file LinearVector.cpp
 * @brief Source of the implementation of the vector linear test equation
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

/// Define small macros allowing to convert to string
#define MAKE_STR_X( _P ) # _P
#define MAKE_STR( _P ) MAKE_STR_X( _P )

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include MAKE_STR( QUICC_MODEL_PATH/Test/LinearVector.hpp )

// Project includes
//
#include "Base/Typedefs.hpp"

namespace QuICC {

namespace Equations {

namespace Test {

   LinearVector::LinearVector(SharedEquationParameters spEqParams)
      : IVectorEquation(spEqParams)
   {
   }

   LinearVector::~LinearVector()
   {
   }

   void LinearVector::setIdentity(const PhysicalNames::Id name)
   {
      // Set the name
      this->setName(name);

      // Set the variable requirements
      this->setRequirements();
   }

   void LinearVector::setCoupling()
   {
      SpectralComponent_range specRange = this->spectralRange();
      SpectralComponent_iterator specIt;
      for(specIt = specRange.first; specIt != specRange.second; ++specIt)
      {
         this->defineCoupling(*specIt, CouplingInformation::PROGNOSTIC, 1, false, false);
      }
   }

   void LinearVector::setRequirements()
   {
      // Add requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(this->name(), FieldRequirement(false, true, true, true));
   }

}
}
}
