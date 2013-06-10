/** \file IVectorDEquation.cpp
 *  \brief Source of the base implementation of a vector diagnostic equation
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "Equations/IVectorDEquation.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Equations {

   IVectorDEquation::IVectorDEquation(SharedEquationParameters spEqParams)
      : IDiagnosticEquation(spEqParams)
   {
   }

   IVectorDEquation::~IVectorDEquation()
   {
   }

   void IVectorDEquation::setUnknown(Datatypes::SharedVectorVariableType spUnknown)
   {
      this->mspUnknown = spUnknown;
   }

   const Datatypes::VectorVariableType& IVectorDEquation::unknown() const
   {
      // Safety assert
      assert(this->mspUnknown);

      return *this->mspUnknown;
   }

   Datatypes::VectorVariableType& IVectorDEquation::rUnknown()
   {
      // Safety assert
      assert(this->mspUnknown);

      return *this->mspUnknown;
   }
}
}
