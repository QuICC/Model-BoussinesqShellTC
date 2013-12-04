/** 
 * @file TauAnnulusChebyshev.hpp
 * @brief Implementation of the Tau conditions for chebyshev based operators for an annulus radius
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef TAUANNULUSCHEBYSHEV_HPP
#define TAUANNULUSCHEBYSHEV_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "StaticAsserts/StaticAssert.hpp"
#include "Base/Typedefs.hpp"
#include "BoundaryCondition/BoundaryCondition.hpp"
#include "SpectralOperators/ITauBoundary.hpp"
#include "Exceptions/Exception.hpp"

namespace GeoMHDiSCC {

namespace Spectral {

   /**
    * @brief Implementation of the Tau conditions for chebyshev based operators for an annulus radius
    */
   class TauAnnulusChebyshev: public ITauBoundary
   {
      public:
         /**
          * @brief Constructor
          *
          * @param c    Boundary condition prefactor
          * @param nN   Size of the Tau basis
          * @param bcs  Vector of boundary conditions
          * @param nEq  Number of equations (only independent of BC for coupled systems)
          */
         TauAnnulusChebyshev(const MHDFloat c, const int nN, const Boundary::BCVector& bcs, const int nEq);

         /**
          * @brief Empty Destructor
          */
         ~TauAnnulusChebyshev();

         /**
          * @brief Constrain matrix with boundary conditions on block of Kronecker product
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
         template <typename TData> const TData& extend(const TData& spec);

         /**
          * @brief Restrict tau coefficients to galerkin size 
          *
          * @param spec Matrix of spectral coefficients
          */
         template <typename TData> const TData& restrict(const TData& spec);
         
      protected:

      private:
         /**
          * @brief Get value at boundary
          *
          * @param pt   boundary point
          */
         virtual Array value(Boundary::BCPosition pt) const;

         /**
          * @brief Get first derivative at boundary
          *
          * @param pt   boundary point
          */
         virtual Array firstDerivative(Boundary::BCPosition pt) const;

         /**
          * @brief Get second derivative at boundary
          *
          * @param pt   boundary point
          */
         virtual Array secondDerivative(Boundary::BCPosition pt) const;

         /**
          * @brief C factor
          *
          * @param n Order of polynial
          */
         MHDFloat c(const int n) const;

         /**
          * @brief 1/C factor
          *
          * @param n Order of polynial
          */
         MHDFloat c_1(const int n) const;

         /**
          * @brief Get array of unit values for left boundary
          */
         Array leftUnit() const;

         /**
          * @brief Get array of unit values for right boundary
          */
         Array rightUnit() const;

   };

   template <typename TData> inline const TData& TauAnnulusChebyshev::extend(const TData& spec)
   {
      return spec;
   }

   template <typename TData> inline const TData& TauAnnulusChebyshev::restrict(const TData& spec)
   {
      return spec;
   }

   template <typename TData> void TauAnnulusChebyshev::constrainKronBlock(TData& mat)
   {
   }

   template <typename TData> bool TauAnnulusChebyshev::constrainKronProduct(TData& mat)
   {
      Debug::StaticAssert<false>();

      return true;
   }

   template <> inline bool TauAnnulusChebyshev::constrainKronProduct<SparseMatrixZ>(SparseMatrixZ& mat)
   {
      if(this->mIsComplex)
      {
         mat =  this->mZTau;
      } else
      {
         mat = this->mRTau.cast<SparseMatrixZ::Scalar>();
      }

      return true;
   }

   template <> inline bool TauAnnulusChebyshev::constrainKronProduct<DecoupledZSparse>(DecoupledZSparse& mat)
   {
      if(this->mIsComplex)
      {
         mat.real() = this->mZTau.real();
         mat.imag() = this->mZTau.imag();
      } else
      {
         mat.real() = this->mRTau;
         mat.imag().resize(this->mRTau.rows(),this->mRTau.cols());
      }

      return true;
   }

   template <> inline bool TauAnnulusChebyshev::constrainKronProduct<SparseMatrix>(SparseMatrix& mat)
   {
      if(this->mIsComplex)
      {
         throw Exception("Can not apply complex stencil to real data");
      } else
      {
         mat = this->mRTau;
      }

      return true;
   }

}
}

#endif // TAUANNULUSCHEBYSHEV_HPP
