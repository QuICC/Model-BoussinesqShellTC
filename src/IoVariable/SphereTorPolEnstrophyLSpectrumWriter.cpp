/** 
 * @file SphereTorPolEnstrophyLSpectrumWriter.cpp
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
#include "IoVariable/SphereTorPolEnstrophyLSpectrumWriter.hpp"

// Project includes
//
#include "Base/MathConstants.hpp"

namespace QuICC {

namespace IoVariable {

   SphereTorPolEnstrophyLSpectrumWriter::SphereTorPolEnstrophyLSpectrumWriter(const std::string& prefix, const std::string& type)
      : ISphericalTorPolEnstrophyLSpectrumWriter(prefix, type)
   {
   }

   SphereTorPolEnstrophyLSpectrumWriter::~SphereTorPolEnstrophyLSpectrumWriter()
   {
   }

   void SphereTorPolEnstrophyLSpectrumWriter::init()
   {
      // Spherical volume: 4/3*pi*r_o^3
      this->mVolume = (4.0/3.0)*Math::PI;

      #if defined QUICC_SPATIALSCHEME_BLFM || defined QUICC_SPATIALSCHEME_WLFM
         this->mHasMOrdering = true;
      #endif //defined QUICC_SPATIALSCHEME_BLFM || defined QUICC_SPATIALSCHEME_WLFM
      #if defined QUICC_SPATIALSCHEME_BLFL || defined QUICC_SPATIALSCHEME_WLFL
         this->mHasMOrdering = false;
      #endif // defined QUICC_SPATIALSCHEME_BLFL || defined QUICC_SPATIALSCHEME_WLFL

      ISphericalTorPolEnstrophyLSpectrumWriter::init();
   }

}
}
