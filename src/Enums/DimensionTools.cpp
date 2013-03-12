/** \file IdToHuman.cpp
 *  \brief Source of enum id to human strings converters
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
