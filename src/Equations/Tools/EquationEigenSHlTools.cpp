/** 
 * @file EquationEigenSHlTools.cpp
 * @brief Source of the tools for schemes with spherical harmonic expansions with l spectral ordering
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Configuration includes
//

// System includes
//

// External includes
//
#include <Eigen/Sparse>

// Class include
//
#include "Equations/Tools/EquationEigenSHlTools.hpp"

// Project includes
//
#include "Base/MathConstants.hpp"

namespace GeoMHDiSCC {

namespace Equations {

namespace EigenSHl {

   int fieldCouplingNMat(const SharedResolution spRes)
   {
      return spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>();
   }

}
}
}
