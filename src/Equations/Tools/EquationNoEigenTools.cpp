/** 
 * @file EquationNoEigenTools.cpp
 * @brief Source of the tools for schemes with no eigen direction
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Configuration includes
//

// System includes
//
#include <limits>

// External includes
//
#include <Eigen/Sparse>

// Class include
//
#include "Equations/Tools/EquationNoEigenTools.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Equations {

namespace NoEigen {

   int fieldCouplingNMat(const SharedResolution spRes)
   {
      return 1;
   }
}
}
}
