/** 
 * @file TauChebyshev.hpp
 * @brief Implementation of the Tau conditions for chebyshev based operators
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef TAUCHEBYSHEV_HPP
#define TAUCHEBYSHEV_HPP

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

namespace GeoMHDiSCC {

namespace Spectral {

   /**
    * @brief Implementation of the Tau conditions for chebyshev based operators
    */
   class TauChebyshev: public ITauBoundary
   {
      public:
         /**
          * @brief Constructor
          *
          * @param nN   Size of the Tau basis
          * @param bcs  Vector of boundary conditions
          * @param nEq  Number of equations (only independent of BC for coupled systems)
          */
         TauChebyshev(const int nN, const Boundary::BCVector& bcs, const int nEq);

         /**
          * @brief Empty Destructor
          */
         ~TauChebyshev();

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

   template <typename TData> inline const TData& TauChebyshev::extend(const TData& spec)
   {
      return spec;
   }

   template <typename TData> inline const TData& TauChebyshev::restrict(const TData& spec)
   {
      return spec;
   }

   template <typename TData> TData TauChebyshev::constrain(const TData& mat)
   {
      Debug::StaticAssert<false>();
   }

   template <> inline SparseMatrixZ TauChebyshev::constrain<SparseMatrixZ>(const SparseMatrixZ& mat)
   {
      return mat + this->mZTau;
   }

   template <> inline SparseMatrix TauChebyshev::constrain<SparseMatrix>(const SparseMatrix& mat)
   {
      return mat + this->mRTau;
   }

}
}

#endif // TAUCHEBYSHEV_HPP
