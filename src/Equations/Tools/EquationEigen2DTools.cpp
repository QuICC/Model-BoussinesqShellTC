/** 
 * @file EquationEigen2DTools.cpp
 * @brief Source of the tools for schemes with two eigen direction
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
#include "Equations/Tools/EquationEigen2DTools.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Equations {

namespace Eigen2D {

   int fieldCouplingNMat(const SharedResolution spRes)
   {
      int nMat = 0;

      for(int i = 0; i < spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++i)
      {
         nMat += spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(i);
      }

      return nMat;
   }

}
}
}
