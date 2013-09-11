/** 
 * @file PolynomialTools.cpp
 * @brief Defines some useful constants and tools for polynomial transforms
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 * @version 0.9.0
 * @date 2013-09-11
 */

// System includes
//
#include <cmath>

// External includes
//

// Class include
//
#include "PolynomialTransforms/PolynomialTools.hpp"

// Project includes
//
#include "IoTools/Formatter.hpp"

namespace GeoMHDiSCC {

namespace Transform {

   const MHDFloat PolynomialTools::STD_DEALIASING = 3.0/2.0;

   int PolynomialTools::dealias(const int size)
   {
      return std::ceil(PolynomialTools::STD_DEALIASING*static_cast<MHDFloat>(size));
   }

}
}
