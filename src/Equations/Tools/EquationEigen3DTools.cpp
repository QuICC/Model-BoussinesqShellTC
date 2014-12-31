/** 
 * @file EquationEigen3DTools.cpp
 * @brief Source of the tools for schemes with three eigen direction
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
#include "Equations/Tools/EquationEigen3DTools.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Equations {

namespace Eigen3D {

   int fieldCouplingNMat(const SharedResolution spRes)
   {
      int nMat = 0;

      for(int i = 0; i < spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++i)
      {
         nMat += spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(i)*spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DATB1D>(i,j);
      }

      return nMat;
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
