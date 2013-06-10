/** \file IScalarDEquation.cpp
 *  \brief Source of the base implementation of a scalar diagnostic equation
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "Equations/IScalarDEquation.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Equations {

   IScalarDEquation::IScalarDEquation(SharedEquationParameters spEqParams)
      : IDiagnosticEquation(spEqParams)
   {
   }

   IScalarDEquation::~IScalarDEquation()
   {
   }

   void IScalarDEquation::setUnknown(Datatypes::SharedScalarVariableType spUnknown)
   {
      this->mspUnknown = spUnknown;
   }

   const Datatypes::ScalarVariableType& IScalarDEquation::unknown() const
   {
      // Safety assert
      assert(this->mspUnknown);

      return *this->mspUnknown;
   }

   Datatypes::ScalarVariableType& IScalarDEquation::rUnknown()
   {
      // Safety assert
      assert(this->mspUnknown);

      return *this->mspUnknown;
   }
}
}
