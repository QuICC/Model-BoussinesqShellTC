/** 
 * @file GalerkinAnnulusChebyshev.hpp
 * @brief Implementation of the Galerkin basis for chebyshev based operators for an annulus radius
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef GALERKINANNULUSCHEBYSHEV_HPP
#define GALERKINANNULUSCHEBYSHEV_HPP

// System includes
//
#include <bitset>

// External includes
//

// Project includes
//
#include "StaticAsserts/StaticAssert.hpp"
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "BoundaryCondition/BoundaryCondition.hpp"
#include "Exceptions/Exception.hpp"

namespace GeoMHDiSCC {

namespace Spectral {

   /**
    * @brief Implementation of the Galerkin basis for chebyshev based operators for an annulus radius
    */
   class GalerkinAnnulusChebyshev
   {
      public:
         /**
          * @brief Constructor
          *
          * @param c    Boundary condition prefactor
          * @param nN   Size of the Tau basis
          * @param bcs  Vector of boundary conditions
          * @param nEq  Number of equations to remove (only independent of BC for coupled systems)
          */
         GalerkinAnnulusChebyshev(const MHDFloat c, const int nN, const Boundary::BCVector& bcs, const int nEq);

         /**
          * @brief Empty Destructor
          */
         ~GalerkinAnnulusChebyshev();

         /**
          * @brief Size of the constained operator
          */
         int nN() const;

         /**
          * @brief Number of boundary rows
          */
         int nBc() const;

         /**
          * @brief Constrain matrix with boundary conditions on block of Kronecker Product
          *
          * @param mat  Matrix operator to constrain
          */
         template <typename TData> void constrainKronBlock(TData& mat);

         /**
          * @brief Constrain matrix with boundary conditions on a kronecker product
          *
          * @param mat  Matrix operator to constrain
          */
         template <typename TData> bool constrainKronProduct(TData& mat);

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
         template <typename TData> TData restrict(const TData& spec);
         
      protected:

      private:
         typedef std::bitset<32> FlagType;

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
          * \{
          */
         static const FlagType ZERO_VALUE_LEFT;  
         static const FlagType ZERO_VALUE_RIGHT;  
         static const FlagType ZERO_D1_LEFT;  
         static const FlagType ZERO_D1_RIGHT;  
         static const FlagType ZERO_D2_LEFT;  
         static const FlagType ZERO_D2_RIGHT;  
         ///\}

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
          * @brief Type of boundary condition
          */
         FlagType mBcId;

         /**
          * @brief Real stencil matrix
          */
         SparseMatrix mRStencil;

         /**
          * @brief Complex stencil matrix
          */
         SparseMatrixZ mZStencil;

   };

   template <typename TData> void GalerkinAnnulusChebyshev::constrainKronBlock(TData& mat)
   {
      Debug::StaticAssert<false>();
   }

   template <typename TData> bool GalerkinAnnulusChebyshev::constrainKronProduct(TData& mat)
   {
      return false;
   }

   template <typename TData> TData GalerkinAnnulusChebyshev::extend(const TData& gal)
   {
      Debug::StaticAssert<false>();
   }

   template <> inline void GalerkinAnnulusChebyshev::constrainKronBlock<SparseMatrixZ>(SparseMatrixZ& mat)
   {
      if(this->mIsComplex)
      {
         mat = (mat*this->mZStencil).bottomRows(this->mN - this->mNeq);
      } else
      {
         mat = (mat*this->mRStencil).bottomRows(this->mN - this->mNeq);
      }
   }

   template <> inline void GalerkinAnnulusChebyshev::constrainKronBlock<DecoupledZSparse>(DecoupledZSparse& mat)
   {
      SparseMatrixZ tmpIn = mat.real().cast<MHDComplex>() + Math::cI*mat.imag();
      SparseMatrixZ tmp;
      if(this->mIsComplex)
      {
         tmp = (tmpIn*this->mZStencil).bottomRows(this->mN - this->mNeq);
      } else
      {
         tmp = (tmpIn*this->mRStencil).bottomRows(this->mN - this->mNeq);
      }
      mat.real() = tmp.real();
      mat.imag() = tmp.imag();
   }

   template <> inline void GalerkinAnnulusChebyshev::constrainKronBlock<SparseMatrix>(SparseMatrix& mat)
   {
      if(this->mIsComplex)
      {
         throw Exception("Can not apply complex stencil to real matrix");
      } else
      {
         mat = (mat*this->mRStencil).bottomRows(this->mN - this->mNeq);
      }
   }

   template <> inline Matrix GalerkinAnnulusChebyshev::extend<Matrix>(const Matrix& gal)
   {
      if(this->mIsComplex)
      {
         throw Exception("Can not apply complex stencil to real data");
         return Matrix();
      } else
      {
         return this->mRStencil*gal;
      }
   }

   template <> inline MatrixZ GalerkinAnnulusChebyshev::extend<MatrixZ>(const MatrixZ& gal)
   {
      if(this->mIsComplex)
      {
         return this->mZStencil*gal;
      } else
      {
         return this->mRStencil*gal;
      }
   }

   template <typename TData> inline TData GalerkinAnnulusChebyshev::restrict(const TData& spec)
   {
      TData tmp = spec.bottomRows(spec.rows()-this->mNeq);
      return tmp;
   }

}
}

#endif // GALERKINANNULUSCHEBYSHEV_HPP
