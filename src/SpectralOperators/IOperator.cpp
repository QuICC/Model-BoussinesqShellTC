/** \file IOperator.cpp
 *  \brief Source of the interface for spectral operators with quasi inverses
 */

// System includes
//
#include <assert.h>

// External includes
//

// Class include
//
#include "SpectralOperators/IOperator.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Spectral {

   IOperator::IOperator(const int basisN)
      : mBasisN(basisN)
   {
   }

   IOperator::~IOperator()
   {
   }

   void IOperator::reset(const int basisN)
   {
      // Set dimension
      this->mBasisN = basisN;
   }

   SparseMatrix IOperator::id(const int p)
   {
      // Create storage for the identity
      SparseMatrix idMat(this->basisN(), this->basisN());

      int colStart;
      int colMax;
      if(p >= 0)
      {
         colStart = p;
         colMax = idMat.cols();
      } else
      {
         colStart = 0;
         colMax = idMat.cols()+p;
      }

      // Create left quasi identity
      idMat.reserve(idMat.cols()-std::abs(p));
      for(int j = colStart; j < colMax; ++j)
      {
         // Create column j
         idMat.startVec(j);

         idMat.insertBack(j,j) = 1.0;
      }
      idMat.finalize(); 

      return idMat;
   }

   int IOperator::basisN() const
   {
      return this->mBasisN;
   }

   DecoupledZSparse IOperator::tau(const DecoupledZMatrix& tauLines, const bool atTop) const
   {
      // Create sparse matrix pair for Tau matrix
      DecoupledZSparse tauMat(std::make_pair(SparseMatrix(this->basisN(), this->basisN()),SparseMatrix(this->basisN(), this->basisN())));

      int nBC;
      // Get the number of boundary conditions (inpute should either have zero size, or both matrices with the same size)
      if(lines.first.size() != 0)
      {
         nBC = lines.first.cols();
      } else
      {
         nBC = lines.second.cols();
      }

      int rowStart;
      int rowMax;
      // Get starting row and max number of rows (depends on atTop)
      if(atTop)
      {
         rowStart = 0;
         rowMax = nBC;
      } else
      {
         rowStart = tauMat.first.rows()-nBC;
         rowMax = tauMat.first.rows();
      }

      // Create real tau lines
      if(lines.first.size() != 0)
      {
         tauMat.first.reserve(lines.first.size());
         for(int j = 0; j < tauMat.first.cols(); ++j)
         {
            // Create column j
            tauMat.first.startVec(j);

            // Loop over tau lines
            for(int i = rowStart; i < rowMax; i++)
            {
               if(lines.first(j,i-rowStart) != 0.0)
               {
                  tauMat.first.insertBack(i,j) = lines.first(j,i-rowStart);
               }
            }
         }
      }
      tauMat.first.finalize(); 

      // Create imaginary tau lines
      if(lines.second.size() != 0)
      {
         tauMat.second.reserve(lines.second.size());
         for(int j = 0; j < tauMat.second.cols(); ++j)
         {
            // Create column j
            tauMat.second.startVec(j);

            // Loop over tau lines
            for(int i = rowStart; i < rowMax; i++)
            {
               if(lines.second(j,i-rowStart) != 0.0)
               {
                  tauMat.second.insertBack(i,j) = lines.second(j,i-rowStart);
               }
            }
         }
      }
      tauMat.second.finalize(); 

      return tauMat;
   }

}
}
