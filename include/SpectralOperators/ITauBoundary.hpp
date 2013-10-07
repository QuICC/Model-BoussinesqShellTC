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
          *
          * @param c    Boundary condition prefactor
          * @param nN   Size of the Tau basis
          * @param nEq  Number of equations to remove (only independent of BC for coupled systems)
          */
         ITauBoundary(const MHDFloat c, const int nN, const int nEq);

         /**
          * @brief Empty Destructor
          */
         ~ITauBoundary();

         /**
          * @brief Size of the constained operator
          */
         int nN() const;

         /**
          * @brief Number of boundary rows
          */
         int nBc() const;
         
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
          * @brief Storage for the prefactor
          */
         MHDFloat mC;

         /**
          * @brief Size of the Tau basis
          */
         int mN;

         /**
          * @brief Number of equations
          */
         int mNeq;

         /**
          * @brief Number of boundrary rows
          */
         int mNbc;

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
