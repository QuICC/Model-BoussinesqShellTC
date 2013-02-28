/** \file BoundaryConditions.hpp
 *  \brief Implementation of the spectral boundary conditions
 *
 *  \mhdBug Needs test
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

namespace GeoMHDiSCC {

namespace Spectral {

   /**
    * @brief Implementation of the spectral boundary conditions
    */
   class BoundaryConditions
   {
      public:
         /**
          * @brief List of avaible boundary conditions
          */
         enum Id {
            /// Boundary value
            VALUE, 
            /// First derivative boundary value
            FIRST_DERIVATIVE,
            /// First derivative boundary value
            SECOND_DERIVATIVE,
            /// Beta complex slope boundary value
            BETA_SLOPE};

         /**
          * @brief Convert boundary condition ids into Tau lines
          *
          * @param bcOp Spectral boundary operator
          * @param bcId Map of boundary condition IDs
          */
         static DecoupledZMatrix tauLines(const IBoundary& bcOp, const std::map<BoundaryConditions::Id,IBoundary::Position>& bcId);
         
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
