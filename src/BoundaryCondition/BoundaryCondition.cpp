/** 
 * @file BoundaryCondition.cpp
 * @brief Source of the abstract boundary condition
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// System includes
//

// External includes
//

// Class include
//
#include "BoundaryCondition/BoundaryCondition.hpp"

// Project includes
//

namespace GeoMHDiSCC {

   namespace Boundary {

      BoundaryCondition::BoundaryCondition(const BCType type, const BCPosition position)
         : type(type), position(position)
      {
      }

      BoundaryCondition::~BoundaryCondition()
      {
      }

   }

}
