/** 
 * @file BoundaryCondition.hpp
 * @brief Implementation of a abstract boundary condition
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef BOUNDARYCONDITION_HPP
#define BOUNDARYCONDITION_HPP

// System includes
//
#include <vector>

// External includes
//

// Project includes
//

namespace GeoMHDiSCC {

namespace Boundary {

   /**
    * @brief List of avaible boundary conditions
    */
   enum BCType {
      /// Boundary value
      VALUE, 
      /// First derivative boundary value
      D1,
      /// First derivative boundary value
      D2,
      /// Beta complex slope boundary value
      BETA_SLOPE
   };

   /**
   * @brief List of possible boundaries
   */
   enum BCPosition {
      LEFT,
      RIGHT
   };

   /**
    * @brief Implementation of the spectral boundary conditions
    */
   class BoundaryCondition
   {
      public:
         /**
          * @brief Constructor
          */
         BoundaryCondition(const BCType type, const BCPosition position);

         /**
          * @brief Empty Destructor
          */
         ~BoundaryCondition();

         /**
          * @brief Type of the condition
          */
         BCType type;

         /**
          * @brief Position of the condition
          */
         BCPosition position;
   };

   /// Typedef for a vector of boundary conditions
   typedef std::vector<BoundaryCondition> BCVector;

}
}

#endif // BOUNDARYCONDITIONS_HPP
