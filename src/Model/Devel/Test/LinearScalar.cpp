/** 
 * @file LinearScalar.cpp
 * @brief Source of the implementation of the scalar linear test equation
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
#include MAKE_STR( QUICC_MODEL_PATH/Test/LinearScalar.hpp )

// Project includes
//
#include "Base/Typedefs.hpp"

namespace QuICC {

namespace Equations {

namespace Test {

   LinearScalar::LinearScalar(SharedEquationParameters spEqParams)
      : IScalarEquation(spEqParams)
   {
   }

   LinearScalar::~LinearScalar()
   {
   }

   void LinearScalar::setIdentity(const PhysicalNames::Id name)
   {
      // Set the name
      this->setName(name);

      // Set the variable requirements
      this->setRequirements();
   }

   void LinearScalar::setCoupling()
   {
      this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::PROGNOSTIC, 1, false, false);
   }

   void LinearScalar::setRequirements()
   {
      // Add requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(this->name(), FieldRequirement(true, true, false, false));
   }

}
}
}
