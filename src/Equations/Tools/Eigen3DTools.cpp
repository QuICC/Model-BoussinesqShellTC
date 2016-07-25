/** 
 * @file Eigen3DTools.cpp
 * @brief Source of the tools for schemes with three eigen direction
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "Equations/Tools/Eigen3DTools.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Equations {

   Eigen3DTools::Eigen3DTools()
   {
   }

   Eigen3DTools::~Eigen3DTools()
   {
   }

   std::vector<MHDFloat> Eigen3DTools::identifyEigs(const SharedResolution& spRes, const int matIdx) const
   {
      throw Exception("Not yet implemented!");
      std::vector<MHDFloat> eigs;

      // Fill eigs somehow

      return eigs;
   }

   int Eigen3DTools::computeNMat(const SharedResolution& spRes) const
   {
      int nMat = 0;

      for(int i = 0; i < spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++i)
      {
         nMat += spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(i)*spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DATB1D>(i,j);
      }

      return nMat;
   }

   void Eigen3DTools::interpretTauN(ArrayI& rTauNs, const int tauSize, const SharedResolution& spRes) const
   {
      rTauNs.setConstant(tauSize);
   }

   void Eigen3DTools::interpretGalerkinN(ArrayI& rGalerkinNs, const int galerkinSize, const SharedResolution& spRes) const
   {
      rGalerkinNs.setConstant(galerkinSize);
   }

   void Eigen3DTools::interpretRhsN(ArrayI& rRhsCols, const int rhsSize, const SharedResolution& spRes) const
   {
      rRhsCols.setConstant(rhsSize);
   }

   void Eigen3DTools::interpretSystemN(ArrayI& rSystemNs, const int systemSize, const SharedResolution& spRes) const
   {
      rSystemNs.setConstant(systemSize);
   }

}
}
