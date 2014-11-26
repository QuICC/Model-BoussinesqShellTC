/** 
 * @file KineticEnergyTags.cpp
 * @brief Source of the definitions and names used by the kinetic energy writer
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// System includes
//

// External includes
//

// Class include
//
#include "IoVariable/KineticEnergyTags.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace IoVariable {

   const std::string KineticEnergyTags::HEADER = "Kinetic energy";

   const std::string KineticEnergyTags::VERSION = "1.0";

   const std::string KineticEnergyTags::BASENAME = "kinetic";

   const std::string KineticEnergyTags::EXTENSION = ".dat";
}
}
