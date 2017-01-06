/** 
 * @file EigenSHlmTools.cpp
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
#include "Equations/Tools/EigenSHlmTools.hpp"

// Project includes
//
#include "Base/MathConstants.hpp"

namespace QuICC {

namespace Equations {

   EigenSHlmTools::EigenSHlmTools()
   {
   }

   EigenSHlmTools::~EigenSHlmTools()
   {
   }

   std::vector<MHDFloat> EigenSHlmTools::identifyEigs(const SharedResolution& spRes, const int matIdx) const
   {
      std::vector<MHDFloat> eigs;

      // Get mode indexes
      ArrayI mode = spRes->cpu()->dim(Dimensions::Transform::TRA1D)->mode(matIdx);

      eigs.push_back(static_cast<MHDFloat>(mode(2)));
      eigs.push_back(static_cast<MHDFloat>(mode(3)));

      return eigs;
   }

   int EigenSHlmTools::computeNMat(const SharedResolution& spRes) const
   {
      int nMat = 0;

      for(int i = 0; i < spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++i)
      {
         nMat += spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(i);
      }

      return nMat;
   }

   void EigenSHlmTools::interpretTauN(ArrayI& rTauNs, const SharedResolution& spRes) const
   {
      // Python setup is sufficient
   }

   void EigenSHlmTools::interpretGalerkinN(ArrayI& rGalerkinNs, const SharedResolution& spRes) const
   {
      // Python setup is sufficient
   }

   void EigenSHlmTools::interpretRhsN(ArrayI& rRhsCols, const SharedResolution& spRes) const
   {
      // Python setup is sufficient
   }

   void EigenSHlmTools::interpretSystemN(ArrayI& rSystemNs, const SharedResolution& spRes) const
   {
      // Python setup is sufficient
   }

}
}
