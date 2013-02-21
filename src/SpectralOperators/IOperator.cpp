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

   IOperator::IOperator(const int polyN)
      : mPolyN(polyN)
   {
   }

   IOperator::~IOperator()
   {
   }

   void IOperator::reset(const int polyN)
   {
      // Set dimension
      this->mPolyN = polyN;
   }

   SparseMatrix IOperator::id(const int p)
   {
      // Create storage for the identity
      SparseMatrix idMat(this->polyN(), this->polyN());

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

   int IOperator::polyN() const
   {
      return this->mPolyN;
   }

   DecoupledZSparse IOperator::createSparseTau(const DecoupledZMatrix& lines, const bool atTop) const
   {
      // Create sparse Tau lines
      DecoupledZSparse tau(std::make_pair(SparseMatrix(this->polyN(), this->polyN()),SparseMatrix(this->polyN(), this->polyN())));

      int nBC;
      if(lines.first.size() != 0)
      {
         nBC = lines.first.cols();
      } else
      {
         nBC = lines.second.cols();
      }

      int rowStart;
      int rowMax;
      if(atTop)
      {
         rowStart = 0;
         rowMax = nBC;
      } else
      {
         rowStart = tau.first.rows()-nBC;
         rowMax = tau.first.rows();
      }

      // Create real tau lines
      if(lines.first.size() != 0)
      {
         tau.first.reserve(lines.first.size());
         for(int j = 0; j < tau.first.cols(); ++j)
         {
            // Create column j
            tau.first.startVec(j);

            // Loop over tau lines
            for(int i = rowStart; i < rowMax; i++)
            {
               if(lines.first(j,i-rowStart) != 0.0)
               {
                  tau.first.insertBack(i,j) = lines.first(j,i-rowStart);
               }
            }
         }
      }
      tau.first.finalize(); 

      // Create imaginary tau lines
      if(lines.second.size() != 0)
      {
         tau.second.reserve(lines.second.size());
         for(int j = 0; j < tau.second.cols(); ++j)
         {
            // Create column j
            tau.second.startVec(j);

            // Loop over tau lines
            for(int i = rowStart; i < rowMax; i++)
            {
               if(lines.second(j,i-rowStart) != 0.0)
               {
                  tau.second.insertBack(i,j) = lines.second(j,i-rowStart);
               }
            }
         }
      }
      tau.second.finalize(); 

      return tau;
   }

}
}
