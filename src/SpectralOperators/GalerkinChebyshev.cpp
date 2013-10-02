/** 
 * @file GalerkinChebyshev.cpp
 * @brief Source of the implementation of the spectral Chebyshev Galerkin operators
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// System includes
//
#include <cassert>

// External includes
//
#include <Eigen/Sparse>

// Class include
//
#include "SpectralOperators/GalerkinChebyshev.hpp"

// Project includes
//
#include "Exceptions/Exception.hpp"

#include <iostream>
namespace GeoMHDiSCC {

namespace Spectral {

   GalerkinChebyshev::GalerkinChebyshev(const int nN, const Boundary::BCVector& bcs, const int nEq)
      : mN(nN), mNeq(nEq)
   {
      this->identifyCondition(bcs);

      this->createStencil();
   }

   GalerkinChebyshev::~GalerkinChebyshev()
   {
   }

   void GalerkinChebyshev::identifyCondition(const Boundary::BCVector& bcs)
   {
      this->mBcId = ZERO_VALUE;
   }

   void GalerkinChebyshev::createStencil()
   {
      if(this->mBcId == ZERO_VALUE_LEFT)
      {
         this->zeroValueLeft(mStencil, this->mN);
      } else if(this->mBcId == ZERO_VALUE_RIGHT)
      {
         this->zeroValueRight(mStencil, this->mN);
      } else if(this->mBcId == ZERO_VALUE)
      {
         this->zeroValue(mStencil, this->mN);
      } else if(this->mBcId == ZERO_D1_LEFT)
      {
         this->zeroD1Left(mStencil, this->mN);
      } else if(this->mBcId == ZERO_D1_RIGHT)
      {
         this->zeroD1Right(mStencil, this->mN);
      } else if(this->mBcId == ZERO_D1)
      {
         this->zeroD1(mStencil, this->mN);
      } else if(this->mBcId == ZERO_D2_LEFT)
      {
         this->zeroD2Left(mStencil, this->mN);
      } else if(this->mBcId == ZERO_D2_RIGHT)
      {
         this->zeroD2Right(mStencil, this->mN);
      } else if(this->mBcId == ZERO_D2)
      {
         this->zeroD2(mStencil, this->mN);
      } else if(this->mBcId == ZERO_VALUED1)
      {
         this->zeroVD1(mStencil, this->mN);
      } else if(this->mBcId == ZERO_VALUED2)
      {
         this->zeroVD2(mStencil, this->mN);
      } else if(this->mBcId == ZERO_D1D2)
      {
         this->zeroD1D2(mStencil, this->mN);
      } else
      {
         throw Exception("Stencil has not been implemented!");
      }
   }

   SparseMatrix GalerkinChebyshev::restrictL(const int nEq, const int cols) const
   {
      assert(cols > 0);
      assert(cols-nEq > 0);

      // Create left restriction matrix
      SparseMatrix matL(cols-nEq, cols);

      matL.reserve(cols-nEq);
      for(int j = nEq; j < matL.cols(); ++j)
      {
         // Create column j
         matL.startVec(j);

         matL.insertBack(j-nEq,j) = 1.0;
      }
      matL.finalize(); 

      return matL;
   }

   void GalerkinChebyshev::zeroValueLeft(SparseMatrix& rStencil, const int cols) const
   {
      std::cerr << "---------- THIS VALUE(-1) STENCIL IS UNTESTED ------------" << std::endl;
      assert(cols > 1);

      rStencil.resize(cols,cols-1);
      rStencil.reserve(2*cols);

      // Fill sparse matrix
      for(int j = 0; j < rStencil.cols(); ++j)
      {
         // Create column j
         rStencil.startVec(j);

         // Create diagonal
         rStencil.insertBack(j,j) = 0.5*GalerkinChebyshev::c(j);

         // Create sub diagonal entry for j+1
         if(j < rStencil.rows()-1)
         {
            rStencil.insertBack(j+1,j) = 0.5;
         }
      }
      rStencil.finalize(); 
   }

   void GalerkinChebyshev::zeroValueRight(SparseMatrix& rStencil, const int cols) const
   {
      std::cerr << "---------- THIS VALUE(-1) STENCIL IS UNTESTED ------------" << std::endl;
      assert(cols > 1);

      rStencil.resize(cols,cols-1);
      rStencil.reserve(2*cols);

      // Fill sparse matrix
      for(int j = 0; j < rStencil.cols(); ++j)
      {
         // Create column j
         rStencil.startVec(j);

         // Create diagonal
         rStencil.insertBack(j,j) = -0.5*GalerkinChebyshev::c(j);

         // Create sub diagonal entry for j+1
         if(j < rStencil.rows()-1)
         {
            rStencil.insertBack(j+1,j) = 0.5;
         }
      }
      rStencil.finalize(); 
   }

   void GalerkinChebyshev::zeroValue(SparseMatrix& rStencil, const int cols) const
   {
      assert(cols > 2);

      rStencil.resize(cols,cols-2);
      rStencil.reserve(2*cols);

      // Fill sparse matrix
      for(int j = 0; j < rStencil.cols(); ++j)
      {
         // Create column j
         rStencil.startVec(j);

         // Create diagonal
         rStencil.insertBack(j,j) = -GalerkinChebyshev::c(j);

         // Create sub diagonal entry for j+2
         if(j < rStencil.rows()-2)
         {
            rStencil.insertBack(j+2,j) = 1.0;
         }
      }
      rStencil.finalize(); 
   }

   void GalerkinChebyshev::zeroD1Left(SparseMatrix& rStencil, const int cols) const
   {
      std::cerr << "---------- THIS D1(-1) STENCIL IS UNTESTED ------------" << std::endl;
      assert(cols > 1);

      rStencil.resize(cols,cols-1);
      rStencil.reserve(2*cols);

      // Fill sparse matrix
      for(int j = 0; j < rStencil.cols(); ++j)
      {
         MHDFloat dj = static_cast<MHDFloat>(j);

         // Create column j
         rStencil.startVec(j);

         // Create diagonal
         rStencil.insertBack(j,j) = -((dj+1.0)*(dj+1.0))*GalerkinChebyshev::c(j)/(2.0*dj*dj + 2.0*dj + 1.0);

         // Create sub diagonal entry for j+1
         if(j < rStencil.rows()-1)
         {
            rStencil.insertBack(j+1,j) = (dj*dj)/(2.0*dj*dj + 2.0*dj + 1.0);
         }
      }
      rStencil.finalize(); 
   }

   void GalerkinChebyshev::zeroD1Right(SparseMatrix& rStencil, const int cols) const
   {
      std::cerr << "---------- THIS D1(-1) STENCIL IS UNTESTED ------------" << std::endl;
      assert(cols > 1);

      rStencil.resize(cols,cols-1);
      rStencil.reserve(2*cols);

      // Fill sparse matrix
      for(int j = 0; j < rStencil.cols(); ++j)
      {
         MHDFloat dj = static_cast<MHDFloat>(j);

         // Create column j
         rStencil.startVec(j);

         // Create diagonal
         rStencil.insertBack(j,j) = ((dj+1.0)*(dj+1.0))*GalerkinChebyshev::c(j)/(2.0*dj*dj + 2.0*dj + 1.0);

         // Create sub diagonal entry for j+1
         if(j < rStencil.rows()-1)
         {
            rStencil.insertBack(j+1,j) = (dj*dj)/(2.0*dj*dj + 2.0*dj + 1.0);
         }
      }
      rStencil.finalize(); 
   }

   void GalerkinChebyshev::zeroD1(SparseMatrix& rStencil, const int cols) const
   {
      assert(cols > 2);

      rStencil.resize(cols,cols-2);
      rStencil.reserve(2*cols);

      // Fill sparse matrix
      for(int j = 0; j < rStencil.cols(); ++j)
      {
         MHDFloat dj = static_cast<MHDFloat>(j);

         // Create column j
         rStencil.startVec(j);

         // Create diagonal
         rStencil.insertBack(j,j) = ((dj+2.0)*(dj+2.0))*GalerkinChebyshev::c(j)/(2.0*dj*dj + 4.0*dj + 4.0);

         // Create sub diagonal entry for j+2
         if(j > 0 && j < rStencil.rows()-2)
         {
            rStencil.insertBack(j+2,j) = -(dj*dj)/(2.0*dj*dj + 4.0*dj + 4.0);
         }
      }
      rStencil.finalize(); 
   }

   void GalerkinChebyshev::zeroD2Left(SparseMatrix& rStencil, const int cols) const
   {
      std::cerr << "---------- THIS D2(-1) STENCIL IS UNTESTED ------------" << std::endl;
      assert(cols > 1);

      rStencil.resize(cols,cols-1);
      rStencil.reserve(2*cols);

      // Fill sparse matrix
      for(int j = 0; j < rStencil.cols(); ++j)
      {
         MHDFloat dj = static_cast<MHDFloat>(j);

         // Create column j
         rStencil.startVec(j);

         // Create diagonal
         rStencil.insertBack(j,j) = ((dj+1.0)*(dj+2.0))*GalerkinChebyshev::c(j)/(2.0*(dj*dj + dj + 1.0));

         // Create sub diagonal entry for j+1
         if(j > 1 && j < rStencil.rows()-1)
         {
            rStencil.insertBack(j+1,j) = ((dj-1.0)*dj)/(2.0*(dj*dj + dj + 1.0));
         }
      }
      rStencil.finalize(); 
   }

   void GalerkinChebyshev::zeroD2Right(SparseMatrix& rStencil, const int cols) const
   {
      std::cerr << "---------- THIS D2(1) STENCIL IS UNTESTED ------------" << std::endl;
      assert(cols > 1);

      rStencil.resize(cols,cols-1);
      rStencil.reserve(2*cols);

      // Fill sparse matrix
      for(int j = 0; j < rStencil.cols(); ++j)
      {
         MHDFloat dj = static_cast<MHDFloat>(j);

         // Create column j
         rStencil.startVec(j);

         // Create diagonal
         rStencil.insertBack(j,j) = ((dj+1.0)*(dj+2.0))*GalerkinChebyshev::c(j)/(2.0*(dj*dj + dj + 1.0));

         // Create sub diagonal entry for j+1
         if(j > 1 && j < rStencil.rows()-1)
         {
            rStencil.insertBack(j+1,j) = -((dj-1.0)*dj)/(2.0*(dj*dj + dj + 1.0));
         }
      }
      rStencil.finalize(); 
   }

   void GalerkinChebyshev::zeroD2(SparseMatrix& rStencil, const int cols) const
   {
      std::cerr << "---------- THIS D2 STENCIL IS UNTESTED ------------" << std::endl;
      assert(cols > 2);

      rStencil.resize(cols,cols-2);
      rStencil.reserve(2*cols);

      // Fill sparse matrix
      for(int j = 0; j < rStencil.cols(); ++j)
      {
         MHDFloat dj = static_cast<MHDFloat>(j);

         // Create column j
         rStencil.startVec(j);

         // Create diagonal
         rStencil.insertBack(j,j) = ((dj+2.0)*(dj+2.0)*(dj+3.0))*GalerkinChebyshev::c(j)/(2.0*(dj+1.0)*(dj*dj + 2.0*dj + 6.0));

         // Create sub diagonal entry for j+2
         if(j > 1 && j < rStencil.rows()-2)
         {
            rStencil.insertBack(j+2,j) = -((dj-1.0)*dj*dj)/(2.0*(dj+1.0)*(dj*dj + 2.0*dj + 6.0));
         }
      }
      rStencil.finalize(); 
   }

   void GalerkinChebyshev::zeroVD1(SparseMatrix& rStencil, const int cols) const
   {
      assert(cols > 4);

      rStencil.resize(cols,cols-4);
      rStencil.reserve(3*cols);

      // Fill sparse matrix
      for(int j = 0; j < rStencil.cols(); ++j)
      {
         MHDFloat dj = static_cast<MHDFloat>(j);

         // Create column j
         rStencil.startVec(j);

         // Create diagonal
         rStencil.insertBack(j,j) = GalerkinChebyshev::c(j);

         // Create sub diagonal entry for j+2
         if(j < rStencil.rows()-2)
         {
            rStencil.insertBack(j+2,j) = -2.0*(dj+2.0)/(dj+3.0);
         }

         // Create sub diagonal entry for j+4
         if(j < rStencil.rows()-4)
         {
            rStencil.insertBack(j+4,j) = (dj+1.0)/(dj+3.0);
         }
      }
      rStencil.finalize(); 
   }

   void GalerkinChebyshev::zeroVD2(SparseMatrix& rStencil, const int cols) const
   {
      assert(cols > 4);

      rStencil.resize(cols,cols-4);
      rStencil.reserve(3*cols);

      // Fill sparse matrix
      for(int j = 0; j < rStencil.cols(); ++j)
      {
         MHDFloat dj = static_cast<MHDFloat>(j);

         // Create column j
         rStencil.startVec(j);

         // Create diagonal
         rStencil.insertBack(j,j) = GalerkinChebyshev::c(j);

         // Create sub diagonal entry for j+2
         if(j < rStencil.rows()-2)
         {
            rStencil.insertBack(j+2,j) = -(2.0*(dj+2.0)*(2.0*dj*dj + 8.0*dj + 15.0))/((dj+3)*(2.0*dj*dj + 12.0*dj + 19.0));
         }

         // Create sub diagonal entry for j+4
         if(j < rStencil.rows()-4)
         {
            rStencil.insertBack(j+4,j) = ((dj+1)*(2.0*dj*dj + 4.0*dj + 3.0))/((dj + 3.0)*(2.0*dj*dj + 12.0*dj + 19.0));
         }
      }
      rStencil.finalize(); 
   }

   void GalerkinChebyshev::zeroD1D2(SparseMatrix& rStencil, const int cols) const
   {
      std::cerr << "---------- THIS D1-D2 STENCIL IS UNTESTED ------------" << std::endl;
      assert(cols > 4);

      rStencil.resize(cols,cols-4);
      rStencil.reserve(3*cols);

      // Fill sparse matrix
      for(int j = 0; j < rStencil.cols(); ++j)
      {
         MHDFloat dj = static_cast<MHDFloat>(j);

         // Create column j
         rStencil.startVec(j);

         // Create diagonal
         rStencil.insertBack(j,j) = ((dj+2.0)*(dj+3.0)*(dj+4.0)*(dj+4.0))*GalerkinChebyshev::c(j)/(4.0*(dj*dj*dj*dj + 8.0*dj*dj*dj*dj + 24.0*dj*dj + 32.0*dj + 24.0));

         // Create sub diagonal entry for j+2
         if(j < rStencil.rows()-2)
         {
            rStencil.insertBack(j+2,j) = -(dj*dj*(dj+4.0)*(dj+4.0))/(2.0*(dj*dj*dj*dj + 8.0*dj*dj*dj*dj + 24.0*dj*dj + 32.0*dj + 24.0));
         }

         // Create sub diagonal entry for j+4
         if(j < rStencil.rows()-4)
         {
            rStencil.insertBack(j+4,j) = (dj*dj*(dj+1.0)*(dj+2.0))/(4.0*(dj*dj*dj*dj + 8.0*dj*dj*dj*dj + 24.0*dj*dj + 32.0*dj + 24.0));
         }
      }
      rStencil.finalize(); 
   }

   MHDFloat GalerkinChebyshev::c(const int n) const
   {
      if(n == 0)
      {
         #ifdef GEOMHDISCC_CHEBYSHEV_HAS_C
            return 1.0;
         #else
            return 2.0;
         #endif

      } else if(n < 0)
      {
         return 0.0;

      } else
      {
         return 1.0;
      }
   }

   MHDFloat GalerkinChebyshev::c_1(const int n) const
   {
      if(n == 0)
      {
         #ifdef GEOMHDISCC_CHEBYSHEV_HAS_C
            return 1.0;
         #else
            return 0.5;
         #endif

      } else if(n < 0)
      {
         return 0.0;

      } else
      {
         return 1.0;
      }
   }


}
}
