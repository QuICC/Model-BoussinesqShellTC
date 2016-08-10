/** 
 * @file IEigenTools.cpp
 * @brief Source of the interface to the eigen direction tools
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
#include "Equations/Tools/IEigenTools.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Equations {

   IEigenTools::IEigenTools()
   {
   }

   IEigenTools::~IEigenTools()
   {
   }

   std::vector<MHDFloat> IEigenTools::getEigs(const SharedResolution& spRes, const int matIdx) const
   {
      return this->identifyEigs(spRes, matIdx);
   }

   int IEigenTools::nMat(const SharedResolution& spRes) const
   {
      return this->computeNMat(spRes);
   }

   void IEigenTools::setTauN(ArrayI& rTauNs, const int tauSize, const SharedResolution& spRes) const
   {
      this->interpretTauN(rTauNs, tauSize, spRes);
   }

   void IEigenTools::setGalerkinN(ArrayI& rGalerkinNs, const int galerkinSize, const SharedResolution& spRes) const
   {
      this->interpretGalerkinN(rGalerkinNs, galerkinSize, spRes);
   }

   void IEigenTools::setRhsN(ArrayI& rRhsCols, const int rhsSize, const SharedResolution& spRes) const
   {
      this->interpretRhsN(rRhsCols, rhsSize, spRes);
   }

   void IEigenTools::setSystemN(ArrayI& rSystemNs, const int systemSize, const SharedResolution& spRes) const
   {
      this->interpretSystemN(rSystemNs, systemSize, spRes);
   }

}
}
