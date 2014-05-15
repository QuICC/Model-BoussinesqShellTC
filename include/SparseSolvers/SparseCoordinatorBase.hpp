/** 
 * @file SparseCoordinatorBase.hpp
 * @brief Implementation of the base for a general sparse solver coordinator
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef SPARSECOORDINATORBASE_HPP
#define SPARSECOORDINATORBASE_HPP

// Configuration includes
//
#include "SmartPointers/SharedPtrMacro.h"

// System includes
//

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Equations/IScalarEquation.hpp"
#include "Equations/IVectorEquation.hpp"

#include <iostream>

namespace GeoMHDiSCC {

namespace Solver {

   /**
    * @brief Implementation of the base for a general sparse solver coordinator
    */
   template <template <class,class> class TSolver, bool TIsComplex> class SparseCoordinatorBase;
}
}

// Include complex field spezialization
#include "SparseSolvers/SparseZCoordinatorBase.hpp"

// Include real field spezialization
#include "SparseSolvers/SparseRCoordinatorBase.hpp"

#endif // SPARSECOORDINATORBASE_HPP
