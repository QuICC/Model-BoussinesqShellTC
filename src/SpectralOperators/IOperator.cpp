/** \file IOperator.cpp
 *  \brief Source of the interface for spectral operators with quasi inverses
 */

// System includes
//
#include <cassert>

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

   SparseMatrix IOperator::id(const int p) const
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
      if(tauLines.first.size() != 0)
      {
         nBC = tauLines.first.cols();
      } else
      {
         nBC = tauLines.second.cols();
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
      if(tauLines.first.size() != 0)
      {
         tauMat.first.reserve(tauLines.first.size());
         for(int j = 0; j < tauMat.first.cols(); ++j)
         {
            // Create column j
            tauMat.first.startVec(j);

            // Loop over tau lines
            for(int i = rowStart; i < rowMax; i++)
            {
               if(tauLines.first(j,i-rowStart) != 0.0)
               {
                  tauMat.first.insertBack(i,j) = tauLines.first(j,i-rowStart);
               }
            }
         }
      }
      tauMat.first.finalize(); 

      // Create imaginary tau lines
      if(tauLines.second.size() != 0)
      {
         tauMat.second.reserve(tauLines.second.size());
         for(int j = 0; j < tauMat.second.cols(); ++j)
         {
            // Create column j
            tauMat.second.startVec(j);

            // Loop over tau tauLines
            for(int i = rowStart; i < rowMax; i++)
            {
               if(tauLines.second(j,i-rowStart) != 0.0)
               {
                  tauMat.second.insertBack(i,j) = tauLines.second(j,i-rowStart);
               }
            }
         }
      }
      tauMat.second.finalize(); 

      return tauMat;
   }

}
}
