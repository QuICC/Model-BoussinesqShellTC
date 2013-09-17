/**
 * @file DimensionTools.hpp
 * @brief Definition of some useful tools for the dimensions enums 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef DIMENSIONTOOLS_HPP
#define DIMENSIONTOOLS_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "Enums/Dimensions.hpp"

namespace GeoMHDiSCC {

namespace Dimensions {

      /**
       * @brief Jump to another transform dimension at runtime
       *
       * @param id      Base id
       * @param step    Size of the dimension jump
       */
      Transform::Id  jump(Transform::Id id, int step);
}
}

#endif // DIMENSIONTOOLS_HPP
