/** 
 * @file SphereScalarEnergyWriter.cpp
 * @brief Source of the implementation of the ASCII spherical harmonics energy calculation for scalar field in a sphere
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Configuration includes
//
#include "Framework/FrameworkMacro.h"

// System includes
//
#include <iomanip>

// External includes
//

// Class include
//
#include "IoVariable/SphereScalarEnergyWriter.hpp"

// Project includes
//
#include "Enums/Dimensions.hpp"
#include "Enums/FieldIds.hpp"

namespace QuICC {

namespace IoVariable {

   SphereScalarEnergyWriter::SphereScalarEnergyWriter(const std::string& prefix, const std::string& type)
      : ISphericalScalarEnergyWriter(prefix, type)
   {
   }

   SphereScalarEnergyWriter::~SphereScalarEnergyWriter()
   {
   }

   void SphereScalarEnergyWriter::init()
   {
      // Normalize by sphere volume: 4/3*pi*r_o^3
      this->mVolume = (4.0/3.0)*Math::PI;

      #if defined QUICC_SPATIALSCHEME_BLFM || defined QUICC_SPATIALSCHEME_WLFM
         this->mHasMOrdering = true;
      #endif //defined QUICC_SPATIALSCHEME_BLFM || defined QUICC_SPATIALSCHEME_WLFM
      #if defined QUICC_SPATIALSCHEME_BLFL || defined QUICC_SPATIALSCHEME_WLFL
         this->mHasMOrdering = false;
      #endif // defined QUICC_SPATIALSCHEME_BLFL || defined QUICC_SPATIALSCHEME_WLFL

      ISphericalScalarEnergyWriter::init();
   }

}
}
