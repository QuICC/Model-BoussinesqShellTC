/** \file IPrognosticEquation.cpp
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
#include "Equations/IPrognosticEquation.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Equations {

   IPrognosticEquation::IPrognosticEquation(SharedEquationParameters spEqParams)
      : IEquation(spEqParams)
   {
   }

   IPrognosticEquation::~IPrognosticEquation()
   {
   }
}
}
