/** 
 * @file EquationEigen1DTools.cpp
 * @brief Source of the tools for schemes with a single eigen direction
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
#include "Equations/Tools/EquationEigen1DTools.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Equations {

namespace Eigen1D {

   int fieldCouplingNMat(const SharedResolution spRes)
   {
      return spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>();
   }

   void interpretTauN(ArrayI& rTauNs, const int tauSize, const SharedResolution spRes)
   {
      rTauNs.setConstant(tauSize);
   }

   void interpretGalerkinN(ArrayI& rGalerkinNs, const int galerkinSize, const SharedResolution spRes)
   {
      rGalerkinNs.setConstant(galerkinSize);
   }

   void interpretRhsN(ArrayI& rRhsCols, const int rhsSize, const SharedResolution spRes)
   {
      rRhsCols.setConstant(rhsSize);
   }

   void interpretSystemN(ArrayI& rSystemNs, const int systemSize, const SharedResolution spRes)
   {
      rSystemNs.setConstant(systemSize);
   }
}
}
}
