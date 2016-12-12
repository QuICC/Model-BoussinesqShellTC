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

   void Eigen3DTools::interpretTauN(ArrayI& rTauNs, const SharedResolution& spRes) const
   {
      // Python setup is sufficient
   }

   void Eigen3DTools::interpretGalerkinN(ArrayI& rGalerkinNs, const SharedResolution& spRes) const
   {
      // Python setup is sufficient
   }

   void Eigen3DTools::interpretRhsN(ArrayI& rRhsCols, const SharedResolution& spRes) const
   {
      // Python setup is sufficient
   }

   void Eigen3DTools::interpretSystemN(ArrayI& rSystemNs, const SharedResolution& spRes) const
   {
      // Python setup is sufficient
   }

}
}
