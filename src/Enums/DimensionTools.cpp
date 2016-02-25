/** 
 * @file DimensionTools.cpp
 * @brief Source of utility tools for the dimension IDs
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
#include "Enums/DimensionTools.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Dimensions {

   Transform::Id jump(Transform::Id id, int step)
   {
      return static_cast<Transform::Id>(static_cast<int>(id)+step);
   }

}
}
