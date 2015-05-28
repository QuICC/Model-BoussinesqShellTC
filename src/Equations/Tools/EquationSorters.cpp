/** 
 * @file EquationSorters.cpp
 * @brief Source of the equation sorting functors
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// System includes
//
#include <cassert>
#include <algorithm>

// External includes
//

// Class include
//
#include "Equations/Tools/EquationSorters.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Equations {

namespace Sorters {

   bool EquationType::operator()(SharedIScalarEquation eqA, SharedIScalarEquation eqB)
   {
      return computeEquationType(eqA) < computeEquationType(eqB);
   }

   bool EquationType::operator()(SharedIVectorEquation eqA, SharedIVectorEquation eqB)
   {
      return computeEquationType(eqA) < computeEquationType(eqB);
   }

   int EquationType::computeEquationType(SharedIScalarEquation eqA)
   {
      return static_cast<int>(eqA->couplingInfo(FieldComponents::Spectral::SCALAR).equationType());
   }

   int EquationType::computeEquationType(SharedIVectorEquation eqA)
   { 
      return static_cast<int>(eqA->couplingInfo(FieldComponents::Spectral::ONE).equationType());
   }
}
}
}
