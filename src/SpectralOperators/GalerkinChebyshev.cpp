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

namespace GeoMHDiSCC {

namespace Spectral {

   SparseMatrix GalerkinChebyshev::constrain(const SparseMatrix& mat, const GalerkinCondition::Id bcId, const int nEq)
   {
      SparseMatrix matR;

      if(bcId == GalerkinCondition::ZERO_VALUE_LEFT)
      {
         GalerkinChebyshev::zeroValueLeft(matR, mat.cols());
      } else if(bcId == GalerkinCondition::ZERO_VALUE_RIGHT)
      {
         GalerkinChebyshev::zeroValueRight(matR, mat.cols());
      } else if(bcId == GalerkinCondition::ZERO_VALUE)
      {
         GalerkinChebyshev::zeroValue(matR, mat.cols());
      } else if(bcId == GalerkinCondition::ZERO_D1_LEFT)
      {
         GalerkinChebyshev::zeroD1Left(matR, mat.cols());
      } else if(bcId == GalerkinCondition::ZERO_D1_RIGHT)
      {
         GalerkinChebyshev::zeroD1Right(matR, mat.cols());
      } else if(bcId == GalerkinCondition::ZERO_D1)
      {
         GalerkinChebyshev::zeroD1(matR, mat.cols());
      } else if(bcId == GalerkinCondition::ZERO_D2_LEFT)
      {
         GalerkinChebyshev::zeroD2Left(matR, mat.cols());
      } else if(bcId == GalerkinCondition::ZERO_D2_RIGHT)
      {
         GalerkinChebyshev::zeroD2Right(matR, mat.cols());
      } else if(bcId == GalerkinCondition::ZERO_D2)
      {
         GalerkinChebyshev::zeroD2(matR, mat.cols());
      } else if(bcId == GalerkinCondition::ZERO_VALUED1)
      {
         GalerkinChebyshev::zeroVD1(matR, mat.cols());
      } else if(bcId == GalerkinCondition::ZERO_VALUED2)
      {
         GalerkinChebyshev::zeroVD2(matR, mat.cols());
      } else if(bcId == GalerkinCondition::ZERO_D1D2)
      {
         GalerkinChebyshev::zeroD1D2(matR, mat.cols());
      } else
      {
         throw Exception("Stencil has not been implemented!");
      }

      return GalerkinChebyshev::restrictL(nEq, mat.cols())*mat*matR;
   }

   Matrix GalerkinChebyshev::extend(const Matrix& gal, const GalerkinCondition::Id bcId, const int nEq)
   {
      SparseMatrix matR;

      if(bcId == GalerkinCondition::ZERO_VALUE_LEFT)
      {
         GalerkinChebyshev::zeroValueLeft(matR, gal.rows()+nEq);
      } else if(bcId == GalerkinCondition::ZERO_VALUE_RIGHT)
      {
         GalerkinChebyshev::zeroValueRight(matR, gal.rows()+nEq);
      } else if(bcId == GalerkinCondition::ZERO_VALUE)
      {
         GalerkinChebyshev::zeroValue(matR, gal.rows()+nEq);
      } else if(bcId == GalerkinCondition::ZERO_D1_LEFT)
      {
         GalerkinChebyshev::zeroD1Left(matR, gal.rows()+nEq);
      } else if(bcId == GalerkinCondition::ZERO_D1_RIGHT)
      {
         GalerkinChebyshev::zeroD1Right(matR, gal.rows()+nEq);
      } else if(bcId == GalerkinCondition::ZERO_D1)
      {
         GalerkinChebyshev::zeroD1(matR, gal.rows()+nEq);
      } else if(bcId == GalerkinCondition::ZERO_D2_LEFT)
      {
         GalerkinChebyshev::zeroD2Left(matR, gal.rows()+nEq);
      } else if(bcId == GalerkinCondition::ZERO_D2_RIGHT)
      {
         GalerkinChebyshev::zeroD2Right(matR, gal.rows()+nEq);
      } else if(bcId == GalerkinCondition::ZERO_D2)
      {
         GalerkinChebyshev::zeroD2(matR, gal.rows()+nEq);
      } else if(bcId == GalerkinCondition::ZERO_VALUED1)
      {
         GalerkinChebyshev::zeroVD1(matR, gal.rows()+nEq);
      } else if(bcId == GalerkinCondition::ZERO_VALUED2)
      {
         GalerkinChebyshev::zeroVD2(matR, gal.rows()+nEq);
      } else if(bcId == GalerkinCondition::ZERO_D1D2)
      {
         GalerkinChebyshev::zeroD1D2(matR, gal.rows()+nEq);
      } else
      {
         throw Exception("Stencil has not been implemented!");
      }

      return matR*gal;
   }

   Matrix GalerkinChebyshev::restrict(const Matrix& spec, const GalerkinCondition::Id bcId, const int nEq)
   {
      return spec.topRows(spec.rows()-nEq);
   }

   SparseMatrix GalerkinChebyshev::restrictL(const int nEq, const int cols)
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

   void GalerkinChebyshev::zeroValueLeft(SparseMatrix& rStencil, const int cols)
   {
      assert(cols > 1);

      rStencil.resize(cols,cols-1);
      rStencil.reserve(2*cols);

      // Fill sparse matrix
      for(int j = 0; j < rStencil.cols(); ++j)
      {
         // Create column j
         rStencil.startVec(j);

         // Create diagonal
         rStencil.insertBack(j,j) = 0.5;

         // Create sub diagonal entry for j+1
         if(j < rStencil.rows()-1)
         {
            rStencil.insertBack(j+1,j) = 0.5;
         }
      }
      rStencil.finalize(); 
   }

   void GalerkinChebyshev::zeroValueRight(SparseMatrix& rStencil, const int cols)
   {
      assert(cols > 1);

      rStencil.resize(cols,cols-1);
      rStencil.reserve(2*cols);

      // Fill sparse matrix
      for(int j = 0; j < rStencil.cols(); ++j)
      {
         // Create column j
         rStencil.startVec(j);

         // Create diagonal
         rStencil.insertBack(j,j) = -0.5;

         // Create sub diagonal entry for j+1
         if(j < rStencil.rows()-1)
         {
            rStencil.insertBack(j+1,j) = 0.5;
         }
      }
      rStencil.finalize(); 
   }

   void GalerkinChebyshev::zeroValue(SparseMatrix& rStencil, const int cols)
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
         rStencil.insertBack(j,j) = -1.0*GalerkinChebyshev::c_1(j);

         // Create sub diagonal entry for j+2
         if(j < rStencil.rows()-2)
         {
            rStencil.insertBack(j+2,j) = 1.0*GalerkinChebyshev::c_1(j);
         }
      }
      rStencil.finalize(); 
   }

   void GalerkinChebyshev::zeroD1Left(SparseMatrix& rStencil, const int cols)
   {
      assert(cols > 1);

      throw Exception("Stencil has not been implemented yet!");
      rStencil.resize(cols,cols-1);
   }

   void GalerkinChebyshev::zeroD1Right(SparseMatrix& rStencil, const int cols)
   {
      assert(cols > 1);

      throw Exception("Stencil has not been implemented yet!");
      rStencil.resize(cols,cols-1);
   }

   void GalerkinChebyshev::zeroD1(SparseMatrix& rStencil, const int cols)
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
         rStencil.insertBack(j,j) = ((dj+2.0)*(dj+2.0))/(2.0*dj*dj + 4.0*dj + 4.0);

         // Create sub diagonal entry for j+2
         if(j > 0 && j < rStencil.rows()-2)
         {
            rStencil.insertBack(j+2,j) = -(dj*dj)/(2.0*dj*dj + 4.0*dj + 4.0);
         }
      }
      rStencil.finalize(); 
   }

   void GalerkinChebyshev::zeroD2Left(SparseMatrix& rStencil, const int cols)
   {
      assert(cols > 1);

      throw Exception("Stencil has not been implemented yet!");
      rStencil.resize(cols,cols-1);
   }

   void GalerkinChebyshev::zeroD2Right(SparseMatrix& rStencil, const int cols)
   {
      assert(cols > 1);

      throw Exception("Stencil has not been implemented yet!");
      rStencil.resize(cols,cols-1);
   }

   void GalerkinChebyshev::zeroD2(SparseMatrix& rStencil, const int cols)
   {
      assert(cols > 2);

      throw Exception("Stencil has not been implemented yet!");
      rStencil.resize(cols,cols-2);
   }

   void GalerkinChebyshev::zeroVD1(SparseMatrix& rStencil, const int cols)
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
         rStencil.insertBack(j,j) = 1.0;

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

   void GalerkinChebyshev::zeroVD2(SparseMatrix& rStencil, const int cols)
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
         rStencil.insertBack(j,j) = 1.0;

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

   void GalerkinChebyshev::zeroD1D2(SparseMatrix& rStencil, const int cols)
   {
      assert(cols > 4);

      throw Exception("Stencil has not been implemented yet!");
      rStencil.resize(cols,cols-4);
   }

   MHDFloat GalerkinChebyshev::c(const int n)
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

   MHDFloat GalerkinChebyshev::c_1(const int n)
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
