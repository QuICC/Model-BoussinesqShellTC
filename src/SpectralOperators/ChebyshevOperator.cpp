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
#include "Exception/Exception.hpp"

namespace GeoMHDiSCC {

namespace Spectral {

   ChebyshevOperator::ChebyshevOperator(const int polyN)
      : IOperator(polyN)
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
      SparseMatrix diffMat(this->polyN(), this->polyN());

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
      SparseMatrix tmp(this->polyN() + pq, this->polyN() + pq);

      // Build the inverse
      this->buildInverse(tmp);
      SparseMatrix high = tmp;

      // Compute right derivative order
      for(int i = 1; i < pq; i++)
      {
         high = tmp*high;
      }

      // Create storage for the inverse
      SparseMatrix invMat(this->polyN(), this->polyN());

      // Create left preudo identity to extract rows
      SparseMatrix idL(this->polyN(), this->polyN() + pq);
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
      SparseMatrix idR(this->polyN() + pq, this->polyN());
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

   EPMFloat ChebyshevOperator::c(const int n) const
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
            mat.insertBack(i,j) = static_cast<EPMFloat>(2*j);
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
            mat.insertBack(j-1,j) = (-1.0/static_cast<EPMFloat>(2*(j-1)));
         }

         // Create sub diagonal entry for j+1
         if(j < mat.rows()-1)
         {
            mat.insertBack(j+1,j) = (1.0/static_cast<EPMFloat>(2*(j+1)));
         }
      }
      mat.finalize(); 
   }

   DecoupledZSparse ChebyshevOperator::tau(const std::map<BoundaryConditions::Id,BoundaryConditions::Position>& bcId, const bool atTop)
   {
      // Map iterator
      std::map<BoundaryConditions::Id,BoundaryConditions::Position>::const_iterator mapIt;

      // Count boundary conditions
      int nBCs = 0;
      for(mapIt = bcId.begin(); mapIt != bcId.end(); ++mapIt)
      {
         if(mapIt->second == BoundaryConditions::BOTH)
         {
            nBCs += 2;
         } else
         {
            nBCs += 1;
         }
      }

      // Flags to now if real and imaginary parts are used
      bool hasReal = false;
      bool hasImag = false;

      // Check for compatible sizes
      assert(this->polyN() >= nBCs);

      // Initialise tau lines matrix
      DecoupledZMatrix tauLines(std::make_pair(Matrix(this->polyN(), nBCs),Matrix(this->polyN(), nBCs)));
      tauLines.first.setZero();
      tauLines.second.setZero();

      Array direct(this->polyN());
      for(int i = 0; i < direct.size(); i++)
      {
         direct(i) = (1.0/this->c(i));
      }

      Array alternate(this->polyN());
      for(int i = 0; i < alternate.size(); i++)
      {
         alternate(i) = (1.0/this->c(i))*std::pow(-1.0,i);
      }

      // Storage for boundary values
      Array val(this->polyN());

      // Create boundary values
      int idx = 0;
      for(mapIt = bcId.begin(); mapIt != bcId.end(); ++mapIt)
      {
         switch(mapIt->first)
         {
            case BoundaryConditions::VALUE:
               if(mapIt->second == BoundaryConditions::BOTH || mapIt->second == BoundaryConditions::TOP)
               {
                  tauLines.first.col(idx) = direct;
                  idx++;
               }
               if(mapIt->second == BoundaryConditions::BOTH || mapIt->second == BoundaryConditions::BOTTOM)
               {
                  tauLines.first.col(idx) = alternate;
                  idx++;
               }
               hasReal = true;
               break;
            case BoundaryConditions::FIRST_DERIVATIVE:
               val.setConstant(0.0);
               for(int i = 1; i < val.size(); i++)
               {
                  val(i) = static_cast<EPMFloat>(i*i);
               }
               if(mapIt->second == BoundaryConditions::BOTH || mapIt->second == BoundaryConditions::TOP)
               {
                  tauLines.first.col(idx) = val;
                  idx++;
               }
               if(mapIt->second == BoundaryConditions::BOTH || mapIt->second == BoundaryConditions::BOTTOM)
               {
                  tauLines.first.col(idx) = val.array()*(-1.0*alternate.array());
                  idx++;
               }
               hasReal = true;
               break;
            case BoundaryConditions::SECOND_DERIVATIVE:
               val.setConstant(0.0);
               for(int i = 2; i < val.size(); i++)
               {
                  val(i) = (1.0/3.0)*static_cast<EPMFloat>(std::pow(i,4) - std::pow(i,2));
               }
               if(mapIt->second == BoundaryConditions::BOTH || mapIt->second == BoundaryConditions::TOP)
               {
                  tauLines.first.col(idx) = val;
                  idx++;
               }
               if(mapIt->second == BoundaryConditions::BOTH || mapIt->second == BoundaryConditions::BOTTOM)
               {
                  tauLines.first.col(idx) = val.array()*alternate.array();
                  idx++;
               }
               hasReal = true;
               break;
            case BoundaryConditions::BETA_SLOPE:
               if(mapIt->second == BoundaryConditions::BOTH || mapIt->second == BoundaryConditions::TOP)
               {
                  tauLines.second.col(idx) = static_cast<EPMFloat>(-this->mSpecIdx(this->mIdx))*direct;
                  idx++;
               }
               if(mapIt->second == BoundaryConditions::BOTH || mapIt->second == BoundaryConditions::BOTTOM)
               {
                  tauLines.second.col(idx) = static_cast<EPMFloat>(this->mSpecIdx(this->mIdx))*alternate;
                  idx++;
               }
               hasImag = true;
               break;
            case BoundaryConditions::FIRST_MODE:
               tauLines.first.col(idx).setConstant(0.0);
               tauLines.first.col(idx)(0) = 1.0;
               idx++;
               hasReal = true;
               break;
            default:
               throw Exception("Unknown boundary condition for Chebyshev tau method");
               break;
         }
      }

      // Clear real matrix if not conditions have been set
      if(!hasReal)
      {
         tauLines.first.resize(0,0);
      }

      // Clear imaginary matrix if not conditions have been set
      if(!hasImag)
      {
         tauLines.second.resize(0,0);
      }

      return this->createSparseTau(tauLines, atTop);
   }

}
}
