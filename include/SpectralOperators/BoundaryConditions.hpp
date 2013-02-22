/** \file BoundaryConditions.hpp
 *  \brief Implementation of the spectral boundary conditions
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
#include "Base/General/Typedefs.hpp"
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
         enum Id {VALUE, FIRST_DERIVATIVE, SECOND_DERIVATIVE, BETA_SLOPE}

         /**
          * @brief Convert boundary condition ids into Tau lines
          *
          * @param bcOp Spectral boundary operator
          * @param bcId Map of boundary condition IDs
          */
         static DecoupledZMatrix BoundaryConditions::tauLines(const IBoundary& bcOp, const std::map<BoundaryConditions::Id,IBoundary::Position>& bcId)
         
      protected:

      private:
         /**
          * @brief Constructor
          *
          * @param basisN   Size of the spectral basis
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
