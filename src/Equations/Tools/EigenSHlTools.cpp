/** 
 * @file EigenSHlTools.cpp
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
#include "Equations/Tools/EigenSHlTools.hpp"

// Project includes
//
#include "Base/MathConstants.hpp"

namespace QuICC {

namespace Equations {

   EigenSHlTools::EigenSHlTools()
   {
   }

   EigenSHlTools::~EigenSHlTools()
   {
   }

   std::vector<MHDFloat> EigenSHlTools::identifyEigs(const SharedResolution& spRes, const int matIdx) const
   {
      std::vector<MHDFloat> eigs;

      eigs.push_back(static_cast<MHDFloat>(spRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(matIdx)));

      return eigs;
   }


   int EigenSHlTools::computeNMat(const SharedResolution& spRes) const
   {
      return spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>();
   }

   void EigenSHlTools::interpretTauN(ArrayI& rTauNs, const SharedResolution& spRes) const
   {
      // Python setup is sufficient
   }

   void EigenSHlTools::interpretGalerkinN(ArrayI& rGalerkinNs, const SharedResolution& spRes) const
   {
      // Python setup is sufficient
   }

   void EigenSHlTools::interpretRhsN(ArrayI& rRhsCols, const SharedResolution& spRes) const
   {
      for(int l = 0; l < spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); l++)
      {
         rRhsCols(l) = spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(l);
      }
   }

   void EigenSHlTools::interpretSystemN(ArrayI& rSystemNs, const SharedResolution& spRes) const
   {
      // Python setup is sufficient
   }

}
}
