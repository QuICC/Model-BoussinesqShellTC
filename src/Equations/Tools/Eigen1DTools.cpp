/** 
 * @file Eigen1DTools.cpp
 * @brief Source of the tools for schemes with a single eigen direction
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
#include "Equations/Tools/Eigen1DTools.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Equations {

   Eigen1DTools::Eigen1DTools()
   {
   }

   Eigen1DTools::~Eigen1DTools()
   {
   }

   std::vector<MHDFloat> Eigen1DTools::identifyEigs(const SharedResolution& spRes, const int matIdx) const
   {
      std::vector<MHDFloat> eigs;

      // Get wave number rescale to box size
      eigs.push_back(spRes->sim()->boxScale(Dimensions::Simulation::SIM2D)*static_cast<MHDFloat>(spRes->cpu()->dim(Dimensions::Transform::TRA1D)->template idx<Dimensions::Data::DATND>(matIdx)));

      return eigs;
   }

   int Eigen1DTools::computeNMat(const SharedResolution& spRes) const
   {
      return spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DATND>();
   }

   void Eigen1DTools::interpretTauN(ArrayI& rTauNs, const int tauSize, const SharedResolution& spRes) const
   {
      rTauNs.setConstant(tauSize);
   }

   void Eigen1DTools::interpretGalerkinN(ArrayI& rGalerkinNs, const int galerkinSize, const SharedResolution& spRes) const
   {
      rGalerkinNs.setConstant(galerkinSize);
   }

   void Eigen1DTools::interpretRhsN(ArrayI& rRhsCols, const int rhsSize, const SharedResolution& spRes) const
   {
      rRhsCols.setConstant(rhsSize);
   }

   void Eigen1DTools::interpretSystemN(ArrayI& rSystemNs, const int systemSize, const SharedResolution& spRes) const
   {
      rSystemNs.setConstant(systemSize);
   }
}
}
