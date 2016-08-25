/** 
 * @file Eigen2DTools.cpp
 * @brief Source of the tools for schemes with two eigen direction
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
#include "Equations/Tools/Eigen2DTools.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Equations {

   Eigen2DTools::Eigen2DTools()
   {
   }

   Eigen2DTools::~Eigen2DTools()
   {
   }

   std::vector<MHDFloat> Eigen2DTools::identifyEigs(const SharedResolution& spRes, const int matIdx) const
   {
      std::vector<MHDFloat> eigs;

      // Get mode indexes
      ArrayI mode = spRes->cpu()->dim(Dimensions::Transform::TRA1D)->mode(matIdx);

      int sN = spRes->sim()->dim(Dimensions::Simulation::SIM2D, Dimensions::Space::SPECTRAL);

      // k2D_
      if(mode(2) < sN/2 + (sN % 2))
      {
         eigs.push_back(spRes->sim()->boxScale(Dimensions::Simulation::SIM2D)*static_cast<MHDFloat>(mode(2)));
      } else
      {
         eigs.push_back(spRes->sim()->boxScale(Dimensions::Simulation::SIM2D)*static_cast<MHDFloat>(mode(2) - sN));
      }

      // k3D_
      eigs.push_back(spRes->sim()->boxScale(Dimensions::Simulation::SIM3D)*static_cast<MHDFloat>(mode(3)));
      
      return eigs;
   }

   int Eigen2DTools::computeNMat(const SharedResolution& spRes) const
   {
      int nMat = 0;

      for(int i = 0; i < spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++i)
      {
         nMat += spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(i);
      }

      return nMat;
   }

   void Eigen2DTools::interpretTauN(ArrayI& rTauNs, const int tauSize, const SharedResolution& spRes) const
   {
      rTauNs.setConstant(tauSize);
   }

   void Eigen2DTools::interpretGalerkinN(ArrayI& rGalerkinNs, const int galerkinSize, const SharedResolution& spRes) const
   {
      rGalerkinNs.setConstant(galerkinSize);
   }

   void Eigen2DTools::interpretRhsN(ArrayI& rRhsCols, const int rhsSize, const SharedResolution& spRes) const
   {
      rRhsCols.setConstant(rhsSize);
   }

   void Eigen2DTools::interpretSystemN(ArrayI& rSystemNs, const int systemSize, const SharedResolution& spRes) const
   {
      rSystemNs.setConstant(systemSize);
   }

}
}
