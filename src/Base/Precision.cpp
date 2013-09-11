/** 
 * @file Precision.cpp
 * @brief Implementation of precision related constants and routines
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 * @version 0.9.0
 * @date 2013-09-11
 */

// System includes
//

// External includes
//

// Class include
//
#include "Base/Precision.hpp"

// Project includes
//

namespace GeoMHDiSCC {
   const internal::MHDFloat Precision::PI = precision::acos(MHD_MP(-1.0));

}
