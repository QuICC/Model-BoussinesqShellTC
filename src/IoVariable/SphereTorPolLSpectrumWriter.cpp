/** 
 * @file SphereTorPolLSpectrumWriter.cpp
 * @brief Source of the implementation of the ASCII spherical harmonics energy spectrum calculation for toroidal/poloidal field in a sphere
 */

// Configuration includes
//
#include "Framework/FrameworkMacro.h"

// System includes
//

// External includes
//

// Class include
//
#include "IoVariable/SphereTorPolLSpectrumWriter.hpp"

// Project includes
//
#include "Base/MathConstants.hpp"

namespace QuICC {

namespace IoVariable {

   SphereTorPolLSpectrumWriter::SphereTorPolLSpectrumWriter(const std::string& prefix, const std::string& type)
      : ISphericalTorPolLSpectrumWriter(prefix, type)
   {
   }

   SphereTorPolLSpectrumWriter::~SphereTorPolLSpectrumWriter()
   {
   }

   void SphereTorPolLSpectrumWriter::init()
   {
      // Spherical shell volume: 4/3*pi*r_o^3
      this->mVolume = (4.0/3.0)*Math::PI;

      #if defined QUICC_SPATIALSCHEME_BLFM || defined QUICC_SPATIALSCHEME_WLFM
         this->mHasMOrdering = true;
      #endif //defined QUICC_SPATIALSCHEME_BLFM || defined QUICC_SPATIALSCHEME_WLFM
      #if defined QUICC_SPATIALSCHEME_BLFL || defined QUICC_SPATIALSCHEME_WLFL
         this->mHasMOrdering = false;
      #endif // defined QUICC_SPATIALSCHEME_BLFL || defined QUICC_SPATIALSCHEME_WLFL

      ISphericalTorPolLSpectrumWriter::init();
   }

}
}
