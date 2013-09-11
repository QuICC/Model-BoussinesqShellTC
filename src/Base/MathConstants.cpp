/** 
 * @file MathConstants.cpp
 * @brief Defines some useful math constants
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 * @version 0.9.0
 * @date 2013-09-11
 */

// System includes
//
#include <cmath>
#include <complex>

// External includes
//

// Class include
//
#include "Base/MathConstants.hpp"

// Project includes
//

namespace GeoMHDiSCC {

   const MHDFloat MathConstants::PI = std::acos(-1);

   const MHDComplex MathConstants::cI = MHDComplex(0.0, 1.0);

}
