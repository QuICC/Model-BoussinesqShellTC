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

   void IEigenTools::setTauN(ArrayI& rTauNs, const SharedResolution& spRes) const
   {
      this->interpretTauN(rTauNs, spRes);
   }

   void IEigenTools::setGalerkinN(ArrayI& rGalerkinNs, const SharedResolution& spRes) const
   {
      this->interpretGalerkinN(rGalerkinNs, spRes);
   }

   void IEigenTools::setRhsN(ArrayI& rRhsCols, const SharedResolution& spRes) const
   {
      this->interpretRhsN(rRhsCols, spRes);
   }

   void IEigenTools::setSystemN(ArrayI& rSystemNs, const SharedResolution& spRes) const
   {
      this->interpretSystemN(rSystemNs, spRes);
   }

}
}
