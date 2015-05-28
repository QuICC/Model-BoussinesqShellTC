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
