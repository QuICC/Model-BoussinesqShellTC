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
#include "BoundaryCondition/BoundaryCondition.hpp"

namespace GeoMHDiSCC {

namespace Spectral {

   /**
    * @brief Implementation of the Galerkin basis for chebyshev based operators
    */
   class GalerkinChebyshev
   {
      public:
         /**
          * @brief Constructor
          *
          * @param nN   Size of the Tau basis
          * @param bs   Vector of boundary conditions
          * @param nEq  Number of equations to remove (only independent of BC for coupled systems)
          */
         GalerkinChebyshev(const int nN, const Boundary::BCVector& bcs, const int nEq);

         /**
          * @brief Empty Destructor
          */
         ~GalerkinChebyshev();

         /**
          * @brief Constrain matrix with boundary conditions
          *
          * @param mat  Matrix operator to constrain
          */
         template <typename TData> TData constrain(const TData& mat);

         /**
          * @brief Extend Galerkin basis coefficients to tau expansion 
          *
          * @param spec Matrix of spectral coefficients
          */
         template <typename TData> TData extend(const TData& spec);

         /**
          * @brief Restrict tau coefficients to galerkin size 
          *
          * @param spec Matrix of spectral coefficients
          */
         template <typename TData> const typename Eigen::Ref<TData>& restrict(const TData& spec);
         
      protected:

      private:
         /**
          * @brief Identify the boundary condition
          */
         void identifyCondition(const Boundary::BCVector& bcs);

         /**
          * @brief Create the stencil matrix
          */
         void createStencil();

         /**
          * @brief Left restriction matrix
          *
          * @param nEq  Number of equations to restrict
          * @param cols Number of columns
          */
         SparseMatrix restrictL(const int nEq, const int cols) const;

         /**
          * @brief Stencil matrix for zero value at the left boundary (-1)
          *
          * @param cols Number of columns
          */
         void zeroValueLeft(SparseMatrix& rStencil, const int cols) const;

         /**
          * @brief Stencil matrix for zero value at the right boundary (+1)
          *
          * @param cols Number of columns
          */
         void zeroValueRight(SparseMatrix& rStencil, const int cols) const;

         /**
          * @brief Stencil matrix for zero value at both boundaries (-1, +1)
          *
          * @param cols Number of columns
          */
         void zeroValue(SparseMatrix& rStencil, const int cols) const;

         /**
          * @brief Stencil matrix for zero first derivative at the left boundary (-1)
          *
          * @param cols Number of columns
          */
         void zeroD1Left(SparseMatrix& rStencil, const int cols) const;

         /**
          * @brief Stencil matrix for zero first derivative at the right boundary (+1)
          *
          * @param cols Number of columns
          */
         void zeroD1Right(SparseMatrix& rStencil, const int cols) const;

         /**
          * @brief Stencil matrix for zero first derivative at both boundaries (-1, +1)
          *
          * @param cols Number of columns
          */
         void zeroD1(SparseMatrix& rStencil, const int cols) const;

         /**
          * @brief Stencil matrix for zero second derivative at the left boundary (-1)
          *
          * @param cols Number of columns
          */
         void zeroD2Left(SparseMatrix& rStencil, const int cols) const;

         /**
          * @brief Stencil matrix for zero second derivative at the right boundary (+1)
          *
          * @param cols Number of columns
          */
         void zeroD2Right(SparseMatrix& rStencil, const int cols) const;

         /**
          * @brief Stencil matrix for zero second derivative at both boundaries (-1, +1)
          *
          * @param cols Number of columns
          */
         void zeroD2(SparseMatrix& rStencil, const int cols) const;

         /**
          * @brief Stencil matrix for zero value and zero first derivative at both boundaries (-1, +1)
          *
          * @param cols Number of columns
          */
         void zeroVD1(SparseMatrix& rStencil, const int cols) const;

         /**
          * @brief Stencil matrix for zero value and zero second derivative at both boundaries (-1, +1)
          *
          * @param cols Number of columns
          */
         void zeroVD2(SparseMatrix& rStencil, const int cols) const;

         /**
          * @brief Stencil matrix for zero first derivative and zero second derivative at both boundaries (-1, +1)
          *
          * @param cols Number of columns
          */
         void zeroD1D2(SparseMatrix& rStencil, const int cols) const;

         /**
          * @brief C factor for Chebyshev polynomials
          *
          * @param n Order of polynial
          */
         MHDFloat c(const int n) const;

         /**
          * @brief 1/C factor for Chebyshev polynomials (to avoid 1/0 problems)
          *
          * @param n Order of polynial
          */
         MHDFloat c_1(const int n) const;

         /**
          * @brief Enum for the type of boundary condition
          */
         enum ConditionId {
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

         /**
          * @brief Size of the Tau basis
          */
         int mN;

         /**
          * @brief Number of equations
          */
         int mNeq;

         /**
          * @brief Type of boundary condition
          */
         ConditionId mBcId;

         /**
          * @brief Stencil matrix
          */
         SparseMatrix mStencil;

   };

   template <typename TData> TData GalerkinChebyshev::constrain(const TData& mat)
   {
      return GalerkinChebyshev::restrictL(this->mNeq, this->mN)*mat*this->mStencil;
   }

   template <typename TData> TData GalerkinChebyshev::extend(const TData& gal)
   {
      return this->mStencil*gal;
   }

   template <typename TData> inline const typename Eigen::Ref<TData>& GalerkinChebyshev::restrict(const TData& spec)
   {
      return spec.bottomRows(spec.rows()-this->mNeq);
   }

}
}

#endif // GALERKINCHEBYSHEV_HPP
