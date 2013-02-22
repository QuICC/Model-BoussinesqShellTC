/** \file IOperator.hpp
 *  \brief Interface for a general spectral operator implementation (direct and quasi inverses)
 */

#ifndef IOPERATOR_HPP
#define IOPERATOR_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"

namespace GeoMHDiSCC {

namespace Spectral {

   /**
    * @brief Interface for a general spectral operator implementation (direct and quasi inverses)
    */
   class IOperator
   {
      public:
         /**
          * @brief Constructor
          *
          * @param basisN   Size of the spectral basis
          */
         IOperator(const int basisN);

         /**
          * @brief Empty Destructor
          */
         virtual ~IOperator();

         /**
          * @brief Reset spectral operators and set new dimension
          *
          * @param basisN   New size of the spectral basis
          */
         void reset(const int basisN);

         /**
          * @brief Get the (qasi) identity matrix of order p
          *
          * The 0th order quasi identity is the standard identity matrix. For positiv p the top rows are removed
          * and for negative p the bottom rows are removed
          *
          * @param p Order of the quasi identity
          */
         virtual SparseMatrix id(const int p);

         /**
          * @brief Get the derivative of order p
          *
          * @param nBC  Number of boundary
          * @param p    Order of the derivative
          */
         virtual SparseMatrix diff(const int nBC, const int p) = 0;

         /**
          * @brief Get the quasi inverse derivative for D^-p D^q
          *
          * @param p    Order of the quasi inverse
          * @param q    Order of the derivative
          */
         virtual SparseMatrix qDiff(const int p, const int q) = 0;

         /**
          * @brief Convert tau lines into sparse matrix pair
          *
          * Any unused component (real, imaginary) will have a zero
          * size if it is not required.
          *
          * @param lines   Tau lines. Matrices have zero size or the same
          * @param atTop   Put tau lines at top of matrix?
          */
         DecoupledZSparse  tau(const DecoupledZMatrix& lines, const bool atTop) const;
         
      protected:
         /**
          * @brief Get size of spectral space
          */
         int basisN() const;

      private:
         /**
          * @brief Size of the polynomial basis
          */
         int mBasisN;
   };

}
}

#endif // IOPERATOR_HPP
