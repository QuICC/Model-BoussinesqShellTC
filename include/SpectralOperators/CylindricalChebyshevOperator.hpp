/** 
 * @file CylindricalChebyshevOperator.hpp
 * @brief Implementation of the spectral operators for the chebyshev basis for a cylindrical radius
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 * @version 0.9.0
 * @date 2013-09-11
 */

#ifndef CYLINDRICALCHEBYSHEVOPERATOR_HPP
#define CYLINDRICALCHEBYSHEVOPERATOR_HPP

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
    * \brief Implementation of the spectral operators for the chebyshev basis for a cylindrical radius
    */
   class CylindricalChebyshevOperator: public IOperator
   {
      public:
         /**
          * @brief Constructor
          *
          * @param basisN   Size of the spectral basis
          */
         CylindricalChebyshevOperator(const int basisN);

         /**
          * @brief Empty Destructor
          */
         ~CylindricalChebyshevOperator();

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

#endif // CYLINDRICALCHEBYSHEVOPERATOR_HPP
