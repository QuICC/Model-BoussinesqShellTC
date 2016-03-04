/** 
 * @file ParityTransformTools.cpp
 * @brief Tools for transform aware of parity
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// System includes
//

// External includes
//

// Class include
//
#include "FastTransforms/ParityTransformTools.hpp"

// Project includes
//

#include <iostream>
namespace GeoMHDiSCC {

namespace Transform {

   void ParityTransformTools::extractParityModes(Matrix& rSelected, const Matrix& data, const MatrixI& info, const int rows)
   {
      assert(info.cols() == 2);
      assert(rSelected.cols() >= info.col(1).sum());

      int k = 0;
      for(int i = 0; i < info.rows(); ++i)
      {
         int j0 = info(i,0);
         int cols = info(i, 1);
         rSelected.block(0, k, rows, cols) = data.block(0, j0, rows, cols);
         k += cols; 
      }
   }

   void ParityTransformTools::setParityModes(Matrix& rData, const Matrix& selected, const MatrixI& info, const int rows)
   {
      assert(info.cols() == 2);
      assert(selected.cols() >= info.col(1).sum());

      int k = 0;
      for(int i = 0; i < info.rows(); ++i)
      {
         int j0 = info(i,0);
         int cols = info(i, 1);
         rData.block(0, j0, rows, cols) = selected.block(0, k, rows, cols);
         k += cols; 
      }
   }

   void ParityTransformTools::addParityModes(Matrix& rData, const Matrix& selected, const MatrixI& info, const int rows)
   {
      assert(info.cols() == 2);
      assert(selected.cols() >= info.col(1).sum());

      int k = 0;
      for(int i = 0; i < info.rows(); ++i)
      {
         int j0 = info(i,0);
         int cols = info(i, 1);
         rData.block(0, j0, rows, cols) += selected.block(0, k, rows, cols);
         k += cols; 
      }
   }

   void ParityTransformTools::scaleParityModes(Matrix& rData, const MatrixI& info, const MHDFloat scale, const int rows)
   {
      assert(info.cols() == 2);

      for(int i = 0; i < info.rows(); ++i)
      {
         int j0 = info(i,0);
         int cols = info(i, 1);
         rData.block(0, j0, rows, cols) *= scale;
      }
   }

   void ParityTransformTools::extractParityModes(Matrix& rSelected, const MatrixZ& data, const bool isReal, const MatrixI& info, const int rows)
   {
      assert(info.cols() == 2);
      assert(rSelected.cols() >= info.col(1).sum());

      int k = 0;
      if(isReal)
      {
         for(int i = 0; i < info.rows(); ++i)
         {
            int j0 = info(i,0);
            int cols = info(i, 1);
            rSelected.block(0, k, rows, cols) = data.block(0, j0, rows, cols).real();
            k += cols; 
         }
      } else
      {
         for(int i = 0; i < info.rows(); ++i)
         {
            int j0 = info(i,0);
            int cols = info(i, 1);
            rSelected.block(0, k, rows, cols) = data.block(0, j0, rows, cols).imag();
            k += cols; 
         }
      }
   }

   void ParityTransformTools::setParityModes(MatrixZ& rData, const Matrix& selected, const bool isReal, const MatrixI& info, const int rows)
   {
      assert(info.cols() == 2);
      assert(selected.cols() >= info.col(1).sum());

      int k = 0;
      if(isReal)
      {
         for(int i = 0; i < info.rows(); ++i)
         {
            int j0 = info(i,0);
            int cols = info(i, 1);
            rData.block(0, j0, rows, cols).real() = selected.block(0, k, rows, cols);
            k += cols; 
         }
      } else
      {
         for(int i = 0; i < info.rows(); ++i)
         {
            int j0 = info(i,0);
            int cols = info(i, 1);
            rData.block(0, j0, rows, cols).imag() = selected.block(0, k, rows, cols);
            k += cols; 
         }
      }
   }

   void ParityTransformTools::addParityModes(MatrixZ& rData, const Matrix& selected, const bool isReal, const MatrixI& info, const int rows)
   {
      assert(info.cols() == 2);
      assert(selected.cols() >= info.col(1).sum());

      int k = 0;
      if(isReal)
      {
         for(int i = 0; i < info.rows(); ++i)
         {
            int j0 = info(i,0);
            int cols = info(i, 1);
            rData.block(0, j0, rows, cols).real() += selected.block(0, k, rows, cols);
            k += cols; 
         }
      } else
      {
         for(int i = 0; i < info.rows(); ++i)
         {
            int j0 = info(i,0);
            int cols = info(i, 1);
            rData.block(0, j0, rows, cols).imag() += selected.block(0, k, rows, cols);
            k += cols; 
         }
      }
   }

   void ParityTransformTools::scaleParityModes(MatrixZ& rData, const bool isReal, const MatrixI& info, const MHDFloat scale, const int rows)
   {
      assert(info.cols() == 2);

      if(isReal)
      {
         for(int i = 0; i < info.rows(); ++i)
         {
            int j0 = info(i,0);
            int cols = info(i, 1);
            rData.block(0, j0, rows, cols).real() *= scale;
         }
      } else
      {
         for(int i = 0; i < info.rows(); ++i)
         {
            int j0 = info(i,0);
            int cols = info(i, 1);
            rData.block(0, j0, rows, cols).imag() *= scale;
         }
      }
   }

   void ParityTransformTools::applyOperator(Matrix& rData, const SparseMatrix& op, const MatrixI& info, const MHDFloat scale, const int rows)
   {
      assert(info.cols() == 2);

      int dataRows = rData.rows(); 
      for(int i = 0; i < info.rows(); ++i)
      {
         int j0 = info(i,0);
         int cols = info(i, 1);
         rData.block(0, j0, rows, cols) = scale*op.topRows(rows)*rData.block(0, j0, dataRows, cols);
      }
   }

   void ParityTransformTools::applyOperator(MatrixZ& rData, const bool isReal, const SparseMatrix& op, const MatrixI& info, const MHDFloat scale, const int rows)
   {
      assert(info.cols() == 2);

      int dataRows = rData.rows(); 
      if(isReal)
      {
         for(int i = 0; i < info.rows(); ++i)
         {
            int j0 = info(i,0);
            int cols = info(i, 1);
            rData.block(0, j0, rows, cols).real() = scale*op.topRows(rows)*rData.block(0, j0, dataRows, cols).real();
         }

      } else
      {
         for(int i = 0; i < info.rows(); ++i)
         {
            int j0 = info(i,0);
            int cols = info(i, 1);
            rData.block(0, j0, rows, cols).imag() = scale*op.topRows(rows)*rData.block(0, j0, dataRows, cols).imag();
         }
      }
   }

   void ParityTransformTools::checkRegularity(const Matrix& data, const int rows)
   {
      Array err = Array::Zero(data.cols());
      for(int i = rows-1; i > 0; --i)
      {
         err.transpose() += std::pow(-1,i+1)*data.row(i);
      }
      err *= 2.0;
      std::cerr << "l = 0: " << err(0) - data(0,0) << " l> 0: " << (err.bottomRows(err.rows()-1).transpose() - data.row(0).rightCols(data.cols()-1)).array().abs().maxCoeff() << std::endl;
   }

   void ParityTransformTools::correctRegularity(Matrix& rData, const int rows)
   {
      rData.row(0).setZero();
      for(int i = rows-1; i > 0; --i)
      {
         rData.row(0) += std::pow(-1,i+1)*rData.row(i);
      }
      rData.row(0) *= 2.0;
   }

}
}
