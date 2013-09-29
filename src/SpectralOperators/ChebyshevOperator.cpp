/** 
 * @file ChebyshevOperator.cpp
 * @brief Source of the implementation of the spectral operator for the Chebyshev basis
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "SpectralOperators/ChebyshevOperator.hpp"

// Project includes
//

#include <iostream>
namespace GeoMHDiSCC {

namespace Spectral {

   ChebyshevOperator::ChebyshevOperator(const int basisN)
      : IOperator(basisN)
   {
   }

   ChebyshevOperator::~ChebyshevOperator()
   {
   }

   SparseMatrix ChebyshevOperator::diff(const int nBC, const int p) const
   {
      assert(p > 0);
      assert(nBC >= 0);

      SparseMatrix mat(this->basisN(), this->basisN());

      // Dispatch to correct routine
      if(p == 1)
      {
         this->buildD1(mat, nBC);
      } else if(p == 2)
      {
         this->buildD2(mat, nBC);
      } else if(p == 4)
      {
         this->buildD4(mat, nBC);
      } else
      {
         this->buildDFromD1(mat, p, nBC);
      }

      return mat;
   }

   SparseMatrix ChebyshevOperator::qDiff(const int p, const int q) const
   {
      assert(p > 0);
      assert(q >= 0);
      assert(p >= q);

      SparseMatrix mat(this->basisN(), this->basisN());

      // Dispatch to the most accurate routine
      if(p == q)
      {
         mat = this->id(p);
      } else if(q == 0)
      {
         if(p == 1)
         {
            this->buildQ1(mat);
         } else if(p == 2)
         {
            this->buildQ2(mat);
         } else if(p == 4)
         {
            this->buildQ4(mat);
         } else
         {
            this->buildQFromQ1(mat,p);
         }
      } else
      {
         if(p == 2 && q == 1)
         {
            this->buildQ2D1(mat);
         } else if(p == 4 && q == 2)
         {
            this->buildQ4D2(mat);
         } else
         {
            this->buildQDProduct(mat, p, q);
         }
      }

      return mat;
   }

   MHDFloat ChebyshevOperator::c(const int n) const
   {
      if(n == 0)
      {
         #ifdef GEOMHDISCC_CHEBYSHEV_HAS_C
            return 2.0;
         #else
            return 1.0;
         #endif

      } else if(n < 0)
      {
         return 0.0;

      } else
      {
         return 1.0;
      }
   }

   MHDFloat ChebyshevOperator::c_1(const int n) const
   {
      if(n == 0)
      {
         #ifdef GEOMHDISCC_CHEBYSHEV_HAS_C
            return 0.5;
         #else
            return 1.0;
         #endif

      } else if(n < 0)
      {
         return 0.0;

      } else
      {
         return 1.0;
      }
   }

   void ChebyshevOperator::buildD1(SparseMatrix& mat, const int nBC) const
   {
      assert(nBC >= 0);

      // Reserve the expected number of nonzero elements
      mat.reserve(2*mat.cols()-1);

      // Fill sparse matrix
      int lastRow = std::min(mat.cols()-1,mat.cols() - nBC);
      for(int j = 1; j < mat.cols(); ++j)
      {
         // Create column j
         mat.startVec(j);
         for(int i = (j-1)%2; i < std::min(j,lastRow); i+=2)
         {
            mat.insertBack(i,j) = static_cast<MHDFloat>(2*j)*this->c_1(i);
         }
      }
      mat.finalize(); 
   }

   void ChebyshevOperator::buildD2(SparseMatrix& mat, const int nBC) const
   {
      assert(nBC >= 0);

      // Reserve the expected number of nonzero elements
      mat.reserve(2*mat.cols()-1);

      // Fill sparse matrix
      int lastRow = std::min(mat.cols()-2,mat.cols() - nBC);
      for(int j = 2; j < mat.cols(); ++j)
      {
         // Create column j
         mat.startVec(j);
         for(int i = j%2; i < std::min(j,lastRow); i+=2)
         {
            mat.insertBack(i,j) = static_cast<MHDFloat>(j*(j*j-i*i))*this->c_1(i);
         }
      }
      mat.finalize(); 
   }

   void ChebyshevOperator::buildD4(SparseMatrix& mat, const int nBC) const
   {
      this->buildDFromD1(mat, 4, nBC);
      
      /// \mhdBug This computation currently requires an ugly expression.
      std::cerr << " !!!!! Using DANGEROUS product to compute matrix" << std::endl;
   }

   void ChebyshevOperator::buildQ1(SparseMatrix& mat) const
   {
      // Reserve the expected number of nonzero elements
      mat.reserve(2);

      // Fill sparse matrix
      for(int j = 0; j < mat.cols()-1; ++j)
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
            mat.insertBack(j+1,j) = (1.0/static_cast<MHDFloat>(2*(j+1)))*this->c(j);
         }
      }
      mat.finalize(); 
   }

   void ChebyshevOperator::buildQ2(SparseMatrix& mat) const
   {
      // Reserve the expected number of nonzero elements
      mat.reserve(3);

      // Fill sparse matrix
      for(int j = 0; j < mat.cols()-2; ++j)
      {
         // Create column j
         mat.startVec(j);

         // Create super diagonal entry for j-2
         if(j > 3)
         {
            mat.insertBack(j-2,j) = (1.0/static_cast<MHDFloat>(4*(j-2)*(j-1)));
         }

         // Create diagonal entry
         if(j > 1)
         {
            mat.insertBack(j,j) = (-1.0/static_cast<MHDFloat>(2*(j*j-1)));
         }

         // Create sub diagonal entry for j+2
         if(j < mat.rows()-2)
         {
            mat.insertBack(j+2,j) = (this->c(j)/static_cast<MHDFloat>(4*(j+2)*(j+1)));
         }
      }
      mat.finalize(); 
   }

   void ChebyshevOperator::buildQ3(SparseMatrix& mat) const
   {
      // Reserve the expected number of nonzero elements
      mat.reserve(4);

      // Fill sparse matrix
      for(int j = 0; j < mat.cols()-3; ++j)
      {
         // Create column j
         mat.startVec(j);

         // Create super diagonal entry for j-3
         if(j > 5)
         {
            mat.insertBack(j-3,j) = (-1.0/static_cast<MHDFloat>(8*(j-3)*(j-2)*(j-1)));
         }

         // Create super diagonal entry for j-1
         if(j > 4)
         {
            mat.insertBack(j-1,j) = (3.0/static_cast<MHDFloat>(8*(j-2)*(j-1)*(j+1)));
         }

         // Create sub diagonal entry for j+1
         if(j < mat.rows()-3)
         {
            mat.insertBack(j+1,j) = (-3.0/static_cast<MHDFloat>(8*(j-1)*(j+1)*(j+2)));
         }

         // Create sub diagonal entry for j+3
         if(j < mat.rows()-3)
         {
            mat.insertBack(j+3,j) = (this->c(j)/static_cast<MHDFloat>(8*(j+1)*(j+2)*(j+3)));
         }
      }
      mat.finalize(); 
   }

   void ChebyshevOperator::buildQ4(SparseMatrix& mat) const
   {
      // Reserve the expected number of nonzero elements
      mat.reserve(5);

      // Fill sparse matrix
      for(int j = 0; j < mat.cols()-4; ++j)
      {
         // Create column j
         mat.startVec(j);

         // Create super diagonal entry for j-4
         if(j > 7)
         {
            mat.insertBack(j-4,j) = (1.0/static_cast<MHDFloat>(16*(j-4)*(j-3)*(j-2)*(j-1)));
         }

         // Create super diagonal entry for j-2
         if(j > 5)
         {
            mat.insertBack(j-2,j) = (-1.0/static_cast<MHDFloat>(4*(j-3)*(j-2)*(j-1)*(j+1)));
         }

         // Create diagonal entry
         if(j > 3)
         {
            mat.insertBack(j,j) = (3.0/static_cast<MHDFloat>(8*(j-2)*(j-1)*(j+1)*(j+2)));
         }

         // Create sub diagonal entry for j+2
         if(j > 1 && j < mat.rows()-4)
         {
            mat.insertBack(j+2,j) = (-1.0/static_cast<MHDFloat>(4*(j-1)*(j+1)*(j+2)*(j+3)));
         }

         // Create sub diagonal entry for j+4
         if(j < mat.rows()-4)
         {
            mat.insertBack(j+4,j) = (this->c(j)/static_cast<MHDFloat>(16*(j+1)*(j+2)*(j+3)*(j+4)));
         }
      }
      mat.finalize(); 
   }

   void ChebyshevOperator::buildDFromD1(SparseMatrix& mat, const int q, const int nBC) const
   {
      assert(nBC >= 0);

      this->buildD1(mat,1);
      SparseMatrix d1 = mat;

      for(int i = 1; i < q; ++i)
      {
         mat = d1*mat;
      }

      if(nBC > q)
      {
         mat = this->id(-nBC)*mat;
      }
   }

   void ChebyshevOperator::buildQFromQ1(SparseMatrix& mat, const int p) const
   {
      assert(p > 0);

      // Increase truncation of operator
      SparseMatrix op(this->basisN() + p, this->basisN() + p);
      this->buildQ1(op);
      SparseMatrix highMat = op;

      for(int i = 1; i < p; ++i)
      {
         highMat = op*highMat;
      }

      // Create left pseudo identity to extract rows
      SparseMatrix idL(this->basisN(), this->basisN() + p);
      idL.reserve(idL.rows()-p);
      for(int j = 0; j < idL.rows(); ++j)
      {
         idL.startVec(j);

         if(j >= p)
         {
            idL.insertBack(j,j) = 1.0;
         }
      }
      idL.finalize(); 

      // Create right pseudo identity to extract rows
      SparseMatrix idR(this->basisN() + p, this->basisN());
      idR.reserve(idR.cols()-p);
      for(int j = 0; j < idR.cols()-p; ++j)
      {
         idR.startVec(j);

         idR.insertBack(j,j) = 1.0;
      }
      idR.finalize(); 

      // Extract upper left corner
      mat = idL * highMat * idR;
   }

   void ChebyshevOperator::buildQ2D1(SparseMatrix& mat) const
   {
      // Reserve the expected number of nonzero elements
      mat.reserve(2);

      // Fill sparse matrix
      for(int j = 0; j < mat.cols()-1; ++j)
      {
         // Create column j
         mat.startVec(j);

         // Create super diagonal entry for j-1
         if(j > 2)
         {
            mat.insertBack(j-1,j) = (-1.0/static_cast<MHDFloat>(2*(j-1)));
         }

         // Create sub diagonal entry for j+1
         if(j > 0 && j < mat.rows()-1)
         {
            mat.insertBack(j+1,j) = (1.0/static_cast<MHDFloat>(2*(j+1)))*this->c(j);
         }
      }

      // Fill in correction to Q1
      for(int j = mat.cols()-1; j < mat.cols(); ++j)
      {
         mat.startVec(j);

         // Create super diagonal entry
         if(j > 3)
         {
            mat.insertBack(j-3,j) = (-static_cast<MHDFloat>(j)/static_cast<MHDFloat>(2*(j-2)*(j-3)));
         }

         // Create diagonal entry
         if(j > 1)
         {
            mat.insertBack(j-1,j) = (static_cast<MHDFloat>(j)/static_cast<MHDFloat>(2*(j-2)*(j-1)));
         }
      }
      mat.finalize(); 
   }

   void ChebyshevOperator::buildQ4D2(SparseMatrix& mat) const
   {
      // Reserve the expected number of nonzero elements
      mat.reserve(3);

      // Fill Q2 matrix
      for(int j = 0; j < mat.cols()-2; ++j)
      {
         // Create column j
         mat.startVec(j);

         // Create super diagonal entry for j-2
         if(j > 5)
         {
            mat.insertBack(j-2,j) = (1.0/static_cast<MHDFloat>(4*(j-2)*(j-1)));
         }

         // Create diagonal entry
         if(j > 3)
         {
            mat.insertBack(j,j) = (-1.0/static_cast<MHDFloat>(2*(j*j-1)));
         }

         // Create sub diagonal entry for j+2
         if(j > 1 && j < mat.rows()-2)
         {
            mat.insertBack(j+2,j) = (this->c(j)/static_cast<MHDFloat>(4*(j+2)*(j+1)));
         }
      }

      // Fill in correction to Q2
      for(int j = mat.cols()-2; j < mat.cols(); ++j)
      {
         // Create column j
         mat.startVec(j);

         // (Super diagonal at +4 for col (N-2))
         mat.insertBack(j-6,j) = (-static_cast<MHDFloat>(j*(j-1))/static_cast<MHDFloat>(4*(j-6)*(j-5)*(j-4)*(j-3)));

         // (Diagonal for col (N-2)) + (super diagonal at +4 for col (N-1))
         mat.insertBack(j-4,j) = (static_cast<MHDFloat>(j*(j-1))/static_cast<MHDFloat>(2*(j-5)*(j-4)*(j-3)*(j-3)) + static_cast<MHDFloat>(j)/static_cast<MHDFloat>(2*(j-4)*(j-3)*(j-3)));

         // (Sub diagonal at - 2 for col (N-2)) + (diagonal for col (N-1))
         mat.insertBack(j-2,j) = (-static_cast<MHDFloat>(j*(j-1))/static_cast<MHDFloat>(4*(j-4)*(j-3)*(j-3)*(j-2)) - static_cast<MHDFloat>(j)/static_cast<MHDFloat>((j-3)*(j-1)*(j-3)));

         // (Sub diagonal at - 2 for col (N-1))
         mat.insertBack(j,j) = (1.0/static_cast<MHDFloat>(2*(j-1)*(j-3)));
      }
      mat.finalize(); 
   }

   void ChebyshevOperator::buildQDProduct(SparseMatrix& mat, const int p, const int q) const
   {
      assert(p > 0);
      assert(q >= 0);
      assert(p >= q);

      SparseMatrix qMat(mat.rows(), mat.cols());
      this->buildQFromQ1(qMat, p);
      SparseMatrix dMat(mat.rows(), mat.cols());
      this->buildDFromD1(dMat, q, p);

      mat = qMat*dMat;
   }

}
}
