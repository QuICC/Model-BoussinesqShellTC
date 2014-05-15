/** 
 * @file TimestepCoordinatorBase.hpp
 * @brief Implementation of a general timestep coordinator base
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef TIMESTEPCOORDINATORBASE_HPP
#define TIMESTEPCOORDINATORBASE_HPP

// Configuration includes
//
#include "SmartPointers/SharedPtrMacro.h"

// System includes
//

// External includes
//

// Project includes
//
#include "Base/MathConstants.hpp"
#include "Enums/ModelOperator.hpp"
#include "SparseSolvers/SparseLinearCoordinatorBase.hpp"
#include "Timesteppers/SparseTimestepper.hpp"
#include "Equations/IScalarEquation.hpp"
#include "Equations/IVectorEquation.hpp"

namespace GeoMHDiSCC {

namespace Timestep {

   /**
    * @brief Implementation of general timestepper base
    */
   template <bool TIsComplex> class TimestepCoordinatorBase;
}
}

// Include complex field spezialization
#include "Timesteppers/TimestepZCoordinatorBase.hpp"

// Include real field spezialization
#include "Timesteppers/TimestepRCoordinatorBase.hpp"

#endif // TIMESTEPCOORDINATORBASE_HPP
