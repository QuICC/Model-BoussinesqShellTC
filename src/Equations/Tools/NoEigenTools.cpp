/** 
 * @file NoEigenTools.cpp
 * @brief Source of the tools for schemes with no eigen direction
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
#include "Equations/Tools/NoEigenTools.hpp"

// Project includes
//

namespace QuICC {

namespace Equations {

   NoEigenTools::NoEigenTools()
   {
   }

   NoEigenTools::~NoEigenTools()
   {
   }

   std::vector<MHDFloat> NoEigenTools::identifyEigs(const SharedResolution& spRes, const int matIdx) const
   {
      std::vector<MHDFloat> eigs;

      return eigs;
   }

   int NoEigenTools::computeNMat(const SharedResolution& spRes) const
   {
      return 1;
   }

   void NoEigenTools::interpretTauN(ArrayI& rTauNs, const int tauSize, const SharedResolution& spRes) const
   {
      rTauNs.setConstant(tauSize);
   }

   void NoEigenTools::interpretGalerkinN(ArrayI& rGalerkinNs, const int galerkinSize, const SharedResolution& spRes) const
   {
      rGalerkinNs.setConstant(galerkinSize);
   }

   void NoEigenTools::interpretRhsN(ArrayI& rRhsCols, const int rhsSize, const SharedResolution& spRes) const
   {
      rRhsCols.setConstant(rhsSize);
   }

   void NoEigenTools::interpretSystemN(ArrayI& rSystemNs, const int systemSize, const SharedResolution& spRes) const
   {
      rSystemNs.setConstant(systemSize);
   }
}
}
