/** 
 * @file ChebyshevOperator.hpp
 * @brief Implementation of the spectral operators for the chebyshev basis
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef CHEBYSHEVOPERATOR_HPP
#define CHEBYSHEVOPERATOR_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "SpectralOperators/IOperator.hpp"

namespace GeoMHDiSCC {

namespace Spectral {

   /**
    * @brief Implementation of the spectral operators for the chebyshev basis
    */
   class ChebyshevOperator: public IOperator
   {
      public:
         /**
          * @brief Constructor
          *
          * @param basisN   Size of the spectral basis
          */
         ChebyshevOperator(const int basisN);

         /**
          * @brief Empty Destructor
          */
         ~ChebyshevOperator();

         /**
          * @brief Get the derivative operator of order q
          *
          * @param nBC  Number of boundary
          * @param q    Order of the derivative
          */
         SparseMatrix diff(const int nBC, const int q) const;

         /**
          * @brief Get the quasi inverse derivative for D^-p D^q
          *
          * @param p    Order of the quasi inverse
          * @param q    Order of the derivative
          */
         SparseMatrix qDiff(const int p, const int q) const;
         
      protected:

      private:
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
          * @brief Build \f$D^{1}\f$ operator through recurrence relation
          *
          * @param mat  Output operator matrix
          */
         void buildD1(SparseMatrix& mat, const int nBC) const;

         /**
          * @brief Build \f$D^{2}\f$ operator through recurrence relation
          *
          * @param mat  Output operator matrix
          */
         void buildD2(SparseMatrix& mat, const int nBC) const;

         /**
          * @brief Build \f$D^{4}\f$ operator through recurrence relation
          *
          * @param mat  Output operator matrix
          */
         void buildD4(SparseMatrix& mat, const int nBC) const;

         /**
          * @brief Build \f$D^{q}\f$ operator by building on \f$D^{1}\f$
          *
          * @param mat  Output operator matrix
          * @param q    Order of the derivative
          */
         void buildDFromD1(SparseMatrix& mat, const int q, const int nBC) const;

         /**
          * @brief Build \f$D^{-1}\f$ operator through recurrence relation
          *
          * @param mat  Output operator matrix
          */
         void buildQ1(SparseMatrix& mat) const;

         /**
          * @brief Build \f$D^{-2}\f$ operator through recurrence relation
          *
          * @param mat  Output operator matrix
          */
         void buildQ2(SparseMatrix& mat) const;

         /**
          * @brief Build \f$D^{-3}\f$ operator through recurrence relation
          *
          * @param mat  Output operator matrix
          */
         void buildQ3(SparseMatrix& mat) const;

         /**
          * @brief Build \f$D^{-4}\f$ operator through recurrence relation
          *
          * @param mat  Output operator matrix
          */
         void buildQ4(SparseMatrix& mat) const;

         /**
          * @brief Build \f$D^{-q}\f$ operator by building on\f$D^{-1}\f$ 
          *
          * @param mat  Output operator matrix
          * @param p    Order of the quasi inverse
          */
         void buildQFromQ1(SparseMatrix& mat, const int p) const;

         /**
          * @brief Build \f$D^{-2}D^{1}\f$ operator through recurrence relation
          *
          * @param mat  Output operator matrix
          */
         void buildQ2D1(SparseMatrix& mat) const;

         /**
          * @brief Build \f$D^{-4}D^{2}\f$ operator through recurrence relation
          *
          * @param mat  Output operator matrix
          */
         void buildQ4D2(SparseMatrix& mat) const;

         /**
          * @brief Build \f$D^{-p}D^{q}\f$ operator by product
          *
          * @param mat  Output operator matrix
          */
         void buildQDProduct(SparseMatrix& mat, const int p , const int q) const;
   };

}
}

#endif // CHEBYSHEVOPERATOR_HPP
