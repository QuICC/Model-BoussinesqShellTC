/** \file IDiagnosticEquation.cpp
 *  \brief Source of building block for the implementation of a time dependend evolution equation
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "Equations/IDiagnosticEquation.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Equations {

   IDiagnosticEquation::IDiagnosticEquation(SharedEquationParameters spEqParams)
      : IEquation(spEqParams)
   {
   }

   IDiagnosticEquation::~IDiagnosticEquation()
   {
   }
}
}
