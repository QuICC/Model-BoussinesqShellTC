/** 
 * @file GalerkinSphereChebyshev.cpp
 * @brief Source of the implementation of the spectral Chebyshev Galerkin operators for a sphere radius
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
#include "SpectralOperators/GalerkinSphereChebyshev.hpp"

// Project includes
//
#include "Exceptions/Exception.hpp"

#include <iostream>
namespace GeoMHDiSCC {

namespace Spectral {

   const GalerkinSphereChebyshev::FlagType GalerkinSphereChebyshev::ZERO_VALUE_LEFT= (1 << 0);
   const GalerkinSphereChebyshev::FlagType GalerkinSphereChebyshev::ZERO_VALUE_RIGHT = (1 << 1);
   const GalerkinSphereChebyshev::FlagType GalerkinSphereChebyshev::ZERO_D1_LEFT = (1 << 2);
   const GalerkinSphereChebyshev::FlagType GalerkinSphereChebyshev::ZERO_D1_RIGHT = (1 << 3);
   const GalerkinSphereChebyshev::FlagType GalerkinSphereChebyshev::ZERO_D2_LEFT = (1 << 4);
   const GalerkinSphereChebyshev::FlagType GalerkinSphereChebyshev::ZERO_D2_RIGHT = (1 << 5);

   GalerkinSphereChebyshev::GalerkinSphereChebyshev(const MHDFloat c, const int nN, const Boundary::BCVector& bcs, const int nEq)
      : mN(nN), mNeq(nEq), mIsComplex(false), mRStencil(0,0), mZStencil(0,0)
   {
      if(c != 1.0)
      {
         throw Exception("Galerkin boundary conditions do not (yet) support prefactors");
      }

      this->identifyCondition(bcs);

      this->createStencil();
   }

   GalerkinSphereChebyshev::~GalerkinSphereChebyshev()
   {
   }

   int GalerkinSphereChebyshev::nN() const
   {
      return this->mN - this->mNeq;
   }

   int GalerkinSphereChebyshev::nBc() const
   {
      return 0;
   }

   void GalerkinSphereChebyshev::identifyCondition(const Boundary::BCVector& bcs)
   {
      this->mIsComplex = false;
      this->mBcId = 0;

      for(Boundary::BCVector::const_iterator it = bcs.begin(); it != bcs.end(); ++it)
      {
         if(it->type == Boundary::VALUE)
         {
            if(it->position == Boundary::LEFT)
            {
               this->mBcId |= ZERO_VALUE_LEFT;
            } else
            {
               this->mBcId |= ZERO_VALUE_RIGHT;
            }
         } else if(it->type == Boundary::D1)
         {
            if(it->position == Boundary::LEFT)
            {
               this->mBcId |= ZERO_D1_LEFT;
            } else
            {
               this->mBcId |= ZERO_D1_RIGHT;
            }
         } else if(it->type == Boundary::D2)
         {
            if(it->position == Boundary::LEFT)
            {
               this->mBcId |= ZERO_D2_LEFT;
            } else
            {
               this->mBcId |= ZERO_D2_RIGHT;
            }
         } else
         {
            throw Exception("Stencil has not been implemented for given combination of boundary conditions!");
         }
      }
   }

   void GalerkinSphereChebyshev::createStencil()
   {
      if(this->mBcId == ZERO_VALUE_LEFT)
      {
         this->zeroValueLeft(mRStencil, this->mN);
      } else if(this->mBcId == ZERO_VALUE_RIGHT)
      {
         this->zeroValueRight(mRStencil, this->mN);
      } else if(this->mBcId == ZERO_D1_LEFT)
      {
         this->zeroD1Left(mRStencil, this->mN);
      } else if(this->mBcId == ZERO_D1_RIGHT)
      {
         this->zeroD1Right(mRStencil, this->mN);
      } else if(this->mBcId == ZERO_D2_LEFT)
      {
         this->zeroD2Left(mRStencil, this->mN);
      } else if(this->mBcId == ZERO_D2_RIGHT)
      {
         this->zeroD2Right(mRStencil, this->mN);
      } else if(this->mBcId == (ZERO_VALUE_LEFT | ZERO_VALUE_RIGHT))
      {
         this->zeroValue(mRStencil, this->mN);
      } else if(this->mBcId == (ZERO_D1_LEFT | ZERO_D1_RIGHT))
      {
         this->zeroD1(mRStencil, this->mN);
      } else if(this->mBcId == (ZERO_D2_LEFT | ZERO_D2_RIGHT))
      {
         this->zeroD2(mRStencil, this->mN);
      } else if(this->mBcId == (ZERO_VALUE_LEFT | ZERO_VALUE_RIGHT | ZERO_D1_LEFT | ZERO_D1_RIGHT))
      {
         this->zeroVD1(mRStencil, this->mN);
      } else if(this->mBcId == (ZERO_VALUE_LEFT | ZERO_VALUE_RIGHT | ZERO_D2_LEFT | ZERO_D2_RIGHT))
      {
         this->zeroVD2(mRStencil, this->mN);
      } else if(this->mBcId == (ZERO_D1_LEFT | ZERO_D1_RIGHT | ZERO_D2_LEFT | ZERO_D2_RIGHT))
      {
         this->zeroD1D2(mRStencil, this->mN);
      } else
      {
         throw Exception("Stencil has not been implemented!");
      }
   }

   SparseMatrix GalerkinSphereChebyshev::restrictL(const int nEq, const int cols) const
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

   void GalerkinSphereChebyshev::zeroValueLeft(SparseMatrix& rStencil, const int cols) const
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
         rStencil.insertBack(j,j) = 0.5*this->c(j);

         // Create sub diagonal entry for j+1
         if(j < rStencil.rows()-1)
         {
            rStencil.insertBack(j+1,j) = 0.5;
         }
      }
      rStencil.finalize(); 
   }

   void GalerkinSphereChebyshev::zeroValueRight(SparseMatrix& rStencil, const int cols) const
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
         rStencil.insertBack(j,j) = -0.5*this->c(j);

         // Create sub diagonal entry for j+1
         if(j < rStencil.rows()-1)
         {
            rStencil.insertBack(j+1,j) = 0.5;
         }
      }
      rStencil.finalize(); 
   }

   void GalerkinSphereChebyshev::zeroValue(SparseMatrix& rStencil, const int cols) const
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
         rStencil.insertBack(j,j) = -this->c(j);

         // Create sub diagonal entry for j+2
         if(j < rStencil.rows()-2)
         {
            rStencil.insertBack(j+2,j) = 1.0;
         }
      }
      rStencil.finalize(); 
   }

   void GalerkinSphereChebyshev::zeroD1Left(SparseMatrix& rStencil, const int cols) const
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
         rStencil.insertBack(j,j) = -((dj+1.0)*(dj+1.0))*this->c(j)/(2.0*dj*dj + 2.0*dj + 1.0);

         // Create sub diagonal entry for j+1
         if(j < rStencil.rows()-1)
         {
            rStencil.insertBack(j+1,j) = (dj*dj)/(2.0*dj*dj + 2.0*dj + 1.0);
         }
      }
      rStencil.finalize(); 
   }

   void GalerkinSphereChebyshev::zeroD1Right(SparseMatrix& rStencil, const int cols) const
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
         rStencil.insertBack(j,j) = ((dj+1.0)*(dj+1.0))*this->c(j)/(2.0*dj*dj + 2.0*dj + 1.0);

         // Create sub diagonal entry for j+1
         if(j < rStencil.rows()-1)
         {
            rStencil.insertBack(j+1,j) = (dj*dj)/(2.0*dj*dj + 2.0*dj + 1.0);
         }
      }
      rStencil.finalize(); 
   }

   void GalerkinSphereChebyshev::zeroD1(SparseMatrix& rStencil, const int cols) const
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
         rStencil.insertBack(j,j) = ((dj+2.0)*(dj+2.0))*this->c(j)/(2.0*dj*dj + 4.0*dj + 4.0);

         // Create sub diagonal entry for j+2
         if(j > 0 && j < rStencil.rows()-2)
         {
            rStencil.insertBack(j+2,j) = -(dj*dj)/(2.0*dj*dj + 4.0*dj + 4.0);
         }
      }
      rStencil.finalize(); 
   }

   void GalerkinSphereChebyshev::zeroD2Left(SparseMatrix& rStencil, const int cols) const
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
         rStencil.insertBack(j,j) = ((dj+1.0)*(dj+2.0))*this->c(j)/(2.0*(dj*dj + dj + 1.0));

         // Create sub diagonal entry for j+1
         if(j > 1 && j < rStencil.rows()-1)
         {
            rStencil.insertBack(j+1,j) = ((dj-1.0)*dj)/(2.0*(dj*dj + dj + 1.0));
         }
      }
      rStencil.finalize(); 
   }

   void GalerkinSphereChebyshev::zeroD2Right(SparseMatrix& rStencil, const int cols) const
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
         rStencil.insertBack(j,j) = ((dj+1.0)*(dj+2.0))*this->c(j)/(2.0*(dj*dj + dj + 1.0));

         // Create sub diagonal entry for j+1
         if(j > 1 && j < rStencil.rows()-1)
         {
            rStencil.insertBack(j+1,j) = -((dj-1.0)*dj)/(2.0*(dj*dj + dj + 1.0));
         }
      }
      rStencil.finalize(); 
   }

   void GalerkinSphereChebyshev::zeroD2(SparseMatrix& rStencil, const int cols) const
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
         rStencil.insertBack(j,j) = ((dj+2.0)*(dj+2.0)*(dj+3.0))*this->c(j)/(2.0*(dj+1.0)*(dj*dj + 2.0*dj + 6.0));

         // Create sub diagonal entry for j+2
         if(j > 1 && j < rStencil.rows()-2)
         {
            rStencil.insertBack(j+2,j) = -((dj-1.0)*dj*dj)/(2.0*(dj+1.0)*(dj*dj + 2.0*dj + 6.0));
         }
      }
      rStencil.finalize(); 
   }

   void GalerkinSphereChebyshev::zeroVD1(SparseMatrix& rStencil, const int cols) const
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
         rStencil.insertBack(j,j) = this->c(j);

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

   void GalerkinSphereChebyshev::zeroVD2(SparseMatrix& rStencil, const int cols) const
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
         rStencil.insertBack(j,j) = this->c(j);

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

   void GalerkinSphereChebyshev::zeroD1D2(SparseMatrix& rStencil, const int cols) const
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
         rStencil.insertBack(j,j) = ((dj+2.0)*(dj+3.0)*(dj+4.0)*(dj+4.0))*this->c(j)/(4.0*(dj*dj*dj*dj + 8.0*dj*dj*dj*dj + 24.0*dj*dj + 32.0*dj + 24.0));

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

   MHDFloat GalerkinSphereChebyshev::c(const int n) const
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

   MHDFloat GalerkinSphereChebyshev::c_1(const int n) const
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
