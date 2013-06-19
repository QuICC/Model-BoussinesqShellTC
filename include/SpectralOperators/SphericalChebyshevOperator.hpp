/** \file SphericalChebyshevOperator.hpp
 *  \brief Implementation of the spectral operators for the chebyshev basis for a spherical radius
 */

#ifndef SPHERICALCHEBYSHEVOPERATOR_HPP
#define SPHERICALCHEBYSHEVOPERATOR_HPP

// System includes
//
#include <set>

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "SpectralOperators/IOperator.hpp"

namespace GeoMHDiSCC {

namespace Spectral {

   /**
    * \brief Implementation of the spectral operators for the chebyshev basis for a spherical radius
    */
   class SphericalChebyshevOperator: public IOperator
   {
      public:
         /**
          * @brief Constructor
          *
          * @param basisN   Size of the spectral basis
          */
         SphericalChebyshevOperator(const int basisN);

         /**
          * @brief Empty Destructor
          */
         ~SphericalChebyshevOperator();

         /**
          * @brief Get the derivative operator of order p
          *
          * @param nBC  Number of boundary
          * @param p    Order of the derivative
          */
         SparseMatrix diff(const int nBC, const int p) const;

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
          * @brief C factor
          *
          * @param n Order of polynial
          */
         MHDFloat c(const int n) const;

         /**
          * @brief Build the derivative matrix
          *
          * @param mat Storage to build the matrix into
          */
         void buildDerivative(SparseMatrix& mat) const;

         /**
          * @brief Build the inverse matrix
          *
          * @param mat Storage to build the matrix into
          */
         void buildInverse(SparseMatrix& mat) const;
   };

}
}

#endif // SPHERICALCHEBYSHEVOPERATOR_HPP
