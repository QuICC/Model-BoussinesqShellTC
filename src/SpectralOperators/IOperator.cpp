/** 
 * @file IOperator.cpp
 * @brief Source of the interface for spectral operators with quasi inverses
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
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

   SparseMatrix IOperator::shiftId(const int p) const
   {
      // Create storage for the identity
      SparseMatrix idMat(this->basisN(), this->basisN());

      int colStart;
      int rowStart;
      int colMax;
      if(p >= 0)
      {
         colStart = 0;
         rowStart = p;
         colMax = idMat.cols()-p;
      } else
      {
         colStart = -p;
         rowStart = 0;
         colMax = idMat.cols();
      }

      // Create left quasi identity
      idMat.reserve(idMat.cols()-std::abs(p));
      for(int j = colStart; j < colMax; ++j)
      {
         // Create column j
         idMat.startVec(j);

         idMat.insertBack(j+rowStart,j) = 1.0;
      }
      idMat.finalize(); 

      return idMat;
   }

   int IOperator::basisN() const
   {
      return this->mBasisN;
   }

}
}
