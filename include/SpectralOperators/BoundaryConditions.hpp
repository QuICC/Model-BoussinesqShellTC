/** 
 * @file BoundaryConditions.hpp
 * @brief Implementation of the spectral boundary conditions
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef BOUNDARYCONDITIONS_HPP
#define BOUNDARYCONDITIONS_HPP

// System includes
//
#include <map>

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "SpectralOperators/IBoundary.hpp"
#include "BoundaryCondition/BoundaryCondition.hpp"

namespace GeoMHDiSCC {

namespace Spectral {

   /**
    * @brief Implementation of the spectral boundary conditions
    */
   class BoundaryConditions
   {
      public:
         /**
          * @brief Convert boundary condition ids into Tau lines
          *
          * @param bcOp Spectral boundary operator
          * @param bcId Map of boundary condition IDs
          */
         static DecoupledZMatrix tauLines(const IBoundary& bcOp, const Boundary::BCVector& bcId);

         /**
          * @brief Convert boundary condition ids into Tau matrices
          *
          * @param bcOp Spectral boundary operator
          * @param bcId Map of boundary condition IDs
          */
         static DecoupledZSparse tauMatrix(const IBoundary& bcOp, const Boundary::BCVector& bcId);
         
      protected:

      private:
         /**
          * @brief Constructor
          */
         BoundaryConditions();

         /**
          * @brief Empty Destructor
          */
         ~BoundaryConditions();
   };

}
}

#endif // BOUNDARYCONDITIONS_HPP
