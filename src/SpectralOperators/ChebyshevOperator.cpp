/** \file ChebyshevOperator.cpp
 *  \brief Source of the implementation of the spectral operator for the Chebyshev basis
 */

// System includes
//
#include <assert.h>

// External includes
//

// Class include
//
#include "SpectralOperators/ChebyshevOperator.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Spectral {

   ChebyshevOperator::ChebyshevOperator(const int basisN)
      : IOperator(basisN)
   {
   }

   ChebyshevOperator::~ChebyshevOperator()
   {
   }

   SparseMatrix ChebyshevOperator::diff(const int nBC, const int p)
   {
      // Check that requested order is possible
      assert(p > 0);

      // Create temporary object
      SparseMatrix diffMat(this->basisN(), this->basisN());

      // Build the derivative
      this->buildDerivative(diffMat);

      // Compute right derivative order
      SparseMatrix op = diffMat;
      for(int i = 1; i < p; i++)
      {
         diffMat = op*diffMat;
      }

      if(nBC == 0)
      {
         return diffMat;
      } else
      {
         return this->id(nBC)*diffMat;
      }
   }

   SparseMatrix ChebyshevOperator::qDiff(const int p, const int q)
   {
      // Check that requested order is possible
      assert(p > 0);
      assert(q >= 0);
      assert(p >= q);

      // Get the effective order of the quasi inverse
      int pq = p - q;

      // Create temporary object
      SparseMatrix tmp(this->basisN() + pq, this->basisN() + pq);

      // Build the inverse
      this->buildInverse(tmp);
      SparseMatrix high = tmp;

      // Compute right derivative order
      for(int i = 1; i < pq; i++)
      {
         high = tmp*high;
      }

      // Create storage for the inverse
      SparseMatrix invMat(this->basisN(), this->basisN());

      // Create left preudo identity to extract rows
      SparseMatrix idL(this->basisN(), this->basisN() + pq);
      idL.reserve(idL.rows()-pq);
      for(int j = 0; j < idL.rows(); ++j)
      {
         // Create column j
         idL.startVec(j);

         // Add diagonal
         if(j >= pq)
         {
            idL.insertBack(j,j) = 1.0;
         }
      }
      idL.finalize(); 

      // Create right preudo identity to extract rows
      SparseMatrix idR(this->basisN() + pq, this->basisN());
      idR.reserve(idR.cols()-pq);
      for(int j = 0; j < idR.cols()-pq; ++j)
      {
         // Create column j
         idR.startVec(j);

         idR.insertBack(j,j) = 1.0;
      }
      idR.finalize(); 

      // Extract upper left corner
      invMat = idL * high * idR;

      if(q == 0)
      {
         return invMat;
      } else if(q == p)
      {
         return this->id(p);
      } else
      {
         return this->id(p) * invMat;
      }
   }

   MHDFloat ChebyshevOperator::c(const int n) const
   {
      if(n == 0)
      {
         return 2.0;

      } else if(n < 0)
      {
         return 0.0;

      } else
      {
         return 1.0;
      }
   }

   void ChebyshevOperator::buildDerivative(SparseMatrix& mat) const
   {
      // Reserve the expected number of nonzero elements
      mat.reserve(2*mat.cols()-1);

      // Fill sparse matrix
      for(int j = 1; j < mat.cols(); ++j)
      {
         // Create column j
         mat.startVec(j);

         // Fill column j
         for(int i = (j-1)%2; i < j; i+=2)
         {
            mat.insertBack(i,j) = static_cast<MHDFloat>(2*j);
         }
      }
      mat.finalize(); 
   }

   void ChebyshevOperator::buildInverse(SparseMatrix& mat) const
   {
      // Reserve the expected number of nonzero elements
      mat.reserve(2*mat.cols()-1);

      // Fill sparse matrix
      for(int j = 0; j < mat.cols(); ++j)
      {
         // Create column j
         mat.startVec(j);

         // Create super diagonal entry for j-1
         if(j > 1)
         {
            mat.insertBack(j-1,j) = (-1.0/static_cast<MHDFloat>(2*(j-1)));
         }

         // Create sub diagonal entry for j+1
         if(j < mat.rows()-1)
         {
            mat.insertBack(j+1,j) = (1.0/static_cast<MHDFloat>(2*(j+1)));
         }
      }
      mat.finalize(); 
   }

}
}