/** 
 * @file EquationEigenSHTools.cpp
 * @brief Source of the tools for schemes with spherical harmonic expansions
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
#include "Equations/Tools/EquationEigenSHTools.hpp"

// Project includes
//
#include "Base/MathConstants.hpp"

namespace GeoMHDiSCC {

namespace Equations {

namespace EigenSH {

   int fieldCouplingNMat(const SharedResolution spRes)
   {
      return spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>();
   }

}
}
}
