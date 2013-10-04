/** 
 * @file SShellChebyshevOperator.hpp
 * @brief Implementation of the spectral operators for the chebyshev basis for a spherical shell radius
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef SSHELLCHEBYSHEVOPERATOR_HPP
#define SSHELLCHEBYSHEVOPERATOR_HPP

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
    * @brief Implementation of the spectral operators for the chebyshev basis for a spherical shell radius
    */
   class SShellChebyshevOperator: public IOperator
   {
      public:
         /**
          * @brief Constructor
          *
          * @param basisN   Size of the spectral basis
          */
         SShellChebyshevOperator(const int basisN);

         /**
          * @brief Empty Destructor
          */
         ~SShellChebyshevOperator();

         /**
          * @brief Get the derivative operator of order p
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
          * @brief C factor
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

#endif // SSHELLCHEBYSHEVOPERATOR_HPP
