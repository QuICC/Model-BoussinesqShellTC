/** 
 * @file ITauBoundary.hpp
 * @brief Implementation of the spectral boundary conditions
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef ITAUBOUNDARY_HPP
#define ITAUBOUNDARY_HPP

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
   class ITauBoundary
   {
      public:
         /**
          * @brief Constructor
          */
         ITauBoundary(const int nN, const int nEq);

         /**
          * @brief Empty Destructor
          */
         ~ITauBoundary();
         
      protected:
         /**
          * @brief Convert boundary condition ids into Tau lines
          *
          * @param bcOp Spectral boundary operator
          * @param bcId Map of boundary condition IDs
          */
         DecoupledZMatrix tauLines(const Boundary::BCVector& bcId) const;

         /**
          * @brief Convert boundary condition ids into Tau matrices
          *
          * @param bcOp Spectral boundary operator
          * @param bcId Map of boundary condition IDs
          */
         void createTauMatrix(const Boundary::BCVector& bcId);

         /**
          * @brief Get value at boundary
          *
          * @param pt   boundary point
          */
         virtual Array value(Boundary::BCPosition pt) const = 0;

         /**
          * @brief Get first derivative at boundary
          *
          * @param pt   boundary point
          */
         virtual Array firstDerivative(Boundary::BCPosition pt) const = 0;

         /**
          * @brief Get second derivative at boundary
          *
          * @param pt   boundary point
          */
         virtual Array secondDerivative(Boundary::BCPosition pt) const = 0;

         /**
          * @brief Size of the Tau basis
          */
         int mN;

         /**
          * @brief Number of equations
          */
         int mNeq;

         /**
          * @brief Tau lines matrix is complex?
          */
         bool mIsComplex;

         /**
          * @brief Storage for a real tau lines matrix
          */
         SparseMatrix mRTau;

         /**
          * @brief Storage for a complex tau lines matrix
          */
         SparseMatrixZ mZTau;

      private:
   };

}
}

#endif // ITAUBOUNDARY_HPP
