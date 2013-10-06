/** 
 * @file BoundaryCoordinator.hpp
 * @brief Implementation of a boundary condition coordinator in multiple dimensions
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef BOUNDARYCOORDINATOR_HPP
#define BOUNDARYCOORDINATOR_HPP

// System includes
//
#include <map>

// External includes
//

// Project includes
//

namespace GeoMHDiSCC {

namespace Boundary {

   /**
    * @brief Implementation of a boundary condition coordinator in multiple dimensions
    */
   template <typename TBcIdx, typename TBc1D, typename TBc2D, typename TBc3D> class BoundaryCoordinator
   {
      public:
         /**
          * @brief Constructor
          */
         BoundaryCoordinator();

         /**
          * @brief Empty Destructor
          */
         ~BoundaryCoordinator();

         /**
          * @brief Add boundary condition to first dimension
          *
          * @param idx  Index of the BC
          * @param bc   BC object
          */
         void add1D(const TBcIdx& idx, const TBc1D& bc);

         /**
          * @brief Add boundary condition to second dimension
          *
          * @param idx  Index of the BC
          * @param bc   BC object
          */
         void add2D(const TBcIdx& idx, const TBc2D& bc);

         /**
          * @brief Add boundary condition to third dimension
          *
          * @param idx  Index of the BC
          * @param bc   BC object
          */
         void add3D(const TBcIdx& idx, const TBc3D& bc);

      protected:

      private:
         /**
          * @brief Boundary conditions for the first dimension
          */
         std::map<TBcIdx, TBc1D> m1D;

         /**
          * @brief Boundary conditions for the first dimension
          */
         std::map<TBcIdx, TBc2D> m2D;

         /**
          * @brief Boundary conditions for the first dimension
          */
         std::map<TBcIdx, TBc3D> m3D;
   };

   template <typename TBcIdx, typename TBc1D, typename TBc2D, typename TBc3D> BoundaryCoordinator<TBcIdx,TBc1D,TBc2D,TBc3D>::BoundaryCoordinator()
   {
   }

   template <typename TBcIdx, typename TBc1D, typename TBc2D, typename TBc3D> BoundaryCoordinator<TBcIdx,TBc1D,TBc2D,TBc3D>::~BoundaryCoordinator()
   {
   }

}
}

#endif // BOUNDARYCOORDINATORS_HPP
