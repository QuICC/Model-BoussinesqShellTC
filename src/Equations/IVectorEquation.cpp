/** \file IVectorEquation.cpp
 *  \brief Source of the base implementation of a vector equation
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "Equations/IVectorEquation.hpp"

// Project includes
//

namespace GeoMHDiSCC {

   IVectorEquation::IVectorEquation(SharedEquationParameters spEqParams)
      : IEvolutionEquation(spEqParams)
   {
   }

   IVectorEquation::~IVectorEquation()
   {
   }
}
