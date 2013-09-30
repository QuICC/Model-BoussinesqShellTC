/** 
 * @file GalerkinChebyshev.hpp
 * @brief Implementation of the Galerkin basis for chebyshev based operators
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef GALERKINCHEBYSHEV_HPP
#define GALERKINCHEBYSHEV_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"

namespace GeoMHDiSCC {

namespace Spectral {

   struct GalerkinCondition 
   {
      enum Id {
         ZERO_VALUE_LEFT,
         ZERO_VALUE_RIGHT,
         ZERO_VALUE,
         ZERO_D1_LEFT,
         ZERO_D1_RIGHT,
         ZERO_D1,
         ZERO_D2_LEFT,
         ZERO_D2_RIGHT,
         ZERO_D2,
         ZERO_VALUED1,
         ZERO_VALUED2,
         ZERO_D1D2
      };
   };

   /**
    * @brief Implementation of the Galerkin basis for chebyshev based operators
    */
   class GalerkinChebyshev
   {
      public:
         /**
          * @brief Compute the 2D Laplacian operator for a periodic 2D box
          *
          * @param mat  Matrix operator to constrain
          * @param bcId ID for the boundary condition to impose
          * @param nEq  Number of equations (only independent of BC for coupled systems)
          */
         static SparseMatrix constrain(const SparseMatrix& mat, const GalerkinCondition::Id bcId, const int nEq);

         /**
          * @brief Extend Galerkin basis coefficients to tau expansion 
          *
          * @param spec Matrix of spectral coefficients
          * @param bcId ID for the boundary condition to impose
          * @param nEq  Number of equations (only independent of BC for coupled systems)
          */
         static Matrix extend(const Matrix& spec, const GalerkinCondition::Id bcId, const int nEq);

         /**
          * @brief Restrict tau coefficients to galerkin size 
          *
          * @param spec Matrix of spectral coefficients
          * @param bcId ID for the boundary condition to impose
          * @param nEq  Number of equations (only independent of BC for coupled systems)
          */
         static Matrix restrict(const Matrix& spec, const GalerkinCondition::Id bcId, const int nEq);
         
      protected:
         /**
          * @brief Left restriction matrix
          *
          * @param nEq  Number of equations to restrict
          * @param cols Number of columns
          */
         static SparseMatrix restrictL(const int nEq, const int cols);

         /**
          * @brief Stencil matrix for zero value at the left boundary (-1)
          *
          * @param cols Number of columns
          */
         static void zeroValueLeft(SparseMatrix& rStencil, const int cols);

         /**
          * @brief Stencil matrix for zero value at the right boundary (+1)
          *
          * @param cols Number of columns
          */
         static void zeroValueRight(SparseMatrix& rStencil, const int cols);

         /**
          * @brief Stencil matrix for zero value at both boundaries (-1, +1)
          *
          * @param cols Number of columns
          */
         static void zeroValue(SparseMatrix& rStencil, const int cols);

         /**
          * @brief Stencil matrix for zero first derivative at the left boundary (-1)
          *
          * @param cols Number of columns
          */
         static void zeroD1Left(SparseMatrix& rStencil, const int cols);

         /**
          * @brief Stencil matrix for zero first derivative at the right boundary (+1)
          *
          * @param cols Number of columns
          */
         static void zeroD1Right(SparseMatrix& rStencil, const int cols);

         /**
          * @brief Stencil matrix for zero first derivative at both boundaries (-1, +1)
          *
          * @param cols Number of columns
          */
         static void zeroD1(SparseMatrix& rStencil, const int cols);

         /**
          * @brief Stencil matrix for zero second derivative at the left boundary (-1)
          *
          * @param cols Number of columns
          */
         static void zeroD2Left(SparseMatrix& rStencil, const int cols);

         /**
          * @brief Stencil matrix for zero second derivative at the right boundary (+1)
          *
          * @param cols Number of columns
          */
         static void zeroD2Right(SparseMatrix& rStencil, const int cols);

         /**
          * @brief Stencil matrix for zero second derivative at both boundaries (-1, +1)
          *
          * @param cols Number of columns
          */
         static void zeroD2(SparseMatrix& rStencil, const int cols);

         /**
          * @brief Stencil matrix for zero value and zero first derivative at both boundaries (-1, +1)
          *
          * @param cols Number of columns
          */
         static void zeroVD1(SparseMatrix& rStencil, const int cols);

         /**
          * @brief Stencil matrix for zero value and zero second derivative at both boundaries (-1, +1)
          *
          * @param cols Number of columns
          */
         static void zeroVD2(SparseMatrix& rStencil, const int cols);

         /**
          * @brief Stencil matrix for zero first derivative and zero second derivative at both boundaries (-1, +1)
          *
          * @param cols Number of columns
          */
         static void zeroD1D2(SparseMatrix& rStencil, const int cols);

         /**
          * @brief C factor for Chebyshev polynomials
          *
          * @param n Order of polynial
          */
         static MHDFloat c(const int n);

         /**
          * @brief 1/C factor for Chebyshev polynomials (to avoid 1/0 problems)
          *
          * @param n Order of polynial
          */
         static MHDFloat c_1(const int n);

      private:
         /**
          * @brief Constructor
          */
         GalerkinChebyshev();

         /**
          * @brief Empty Destructor
          */
         ~GalerkinChebyshev();
   };

}
}

#endif // GALERKINCHEBYSHEV_HPP
