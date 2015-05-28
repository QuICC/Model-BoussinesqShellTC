/**
 * @file MatrixOperationsInternal.hpp
 * @brief Useful methods for the DecoupledComplex type
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef MATRIXOPERATIONSINTERNAL_HPP
#define MATRIXOPERATIONSINTERNAL_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Exceptions/Exception.hpp"

namespace GeoMHDiSCC {

namespace Datatypes  {

   namespace internal
   {
      //
      // Matrix products
      //
      
      void setMatrixProduct(Matrix& rField, const int start, const SparseMatrix& mat, const Matrix& rhs);

      void setMatrixProduct(Matrix& rField, const int start, const SparseMatrixZ& mat, const Matrix& rhs);

      void setMatrixProduct(MatrixZ& rField, const int start, const SparseMatrix& mat, const MatrixZ& rhs);

      void setMatrixProduct(MatrixZ& rField, const int start, const SparseMatrixZ& mat, const Eigen::Ref<const MatrixZ>& rhs);

      void setMatrixProduct(DecoupledZMatrix& rField, const int start, const SparseMatrix& mat, const Eigen::Ref<const MatrixZ>& rhs);

      void setMatrixProduct(DecoupledZMatrix& rField, const int start, const SparseMatrix& mat, const Eigen::Ref<const Matrix>& rhsReal, const Eigen::Ref<const Matrix>& rhsImag);

      void setMatrixProduct(DecoupledZMatrix& rField, const int start, const SparseMatrixZ& mat, const Eigen::Ref<const Matrix>& rhsReal, const Eigen::Ref<const Matrix>& rhsImag);

      void setMatrixProduct(DecoupledZMatrix& rField, const int start, const SparseMatrixZ& mat, const Eigen::Ref<const MatrixZ>& rhs);
      
      void addMatrixProduct(Matrix& rField, const int start, const SparseMatrix& mat, const Matrix& rhs);

      void addMatrixProduct(Matrix& rField, const int start, const SparseMatrixZ& mat, const Matrix& rhs);

      void addMatrixProduct(MatrixZ& rField, const int start, const SparseMatrix& mat, const MatrixZ& rhs);

      void addMatrixProduct(MatrixZ& rField, const int start, const SparseMatrixZ& mat, const Eigen::Ref<const MatrixZ>& rhs);

      void addMatrixProduct(DecoupledZMatrix& rField, const int start, const SparseMatrix& mat, const Eigen::Ref<const MatrixZ>& rhs);

      void addMatrixProduct(DecoupledZMatrix& rField, const int start, const SparseMatrix& mat, const Eigen::Ref<const Matrix>& rhsReal, const Eigen::Ref<const Matrix>& rhsImag);

      void addMatrixProduct(DecoupledZMatrix& rField, const int start, const SparseMatrixZ& mat, const Eigen::Ref<const Matrix>& rhsReal, const Eigen::Ref<const Matrix>& rhsImag);

      void addMatrixProduct(DecoupledZMatrix& rField, const int start, const SparseMatrixZ& mat, const Eigen::Ref<const MatrixZ>& rhs);

      //
      // Top block operations
      //

      void setTopBlock(MatrixZ& rField, const int start, const int rows, const MatrixZ& rhs);

      void setTopBlock(Matrix& rField, const int start, const int rows, const Matrix& rhs);

      void setTopBlock(DecoupledZMatrix& rField, const int start, const int rows, const DecoupledZMatrix& rhs);

      /*
       * Inline definitions below
       */

      inline void setMatrixProduct(Matrix& rField, const int start, const SparseMatrixZ& mat, const Matrix& rhs)
      {
         throw Exception("Tried to use complex operator on real fields");
      }

      inline void setMatrixProduct(Matrix& rField, const int start, const SparseMatrix& mat, const Matrix& rhs)
      {
         rField.block(start, 0, mat.rows(), rField.cols()) = mat*rhs;
      }

      inline void setMatrixProduct(MatrixZ& rField, const int start, const SparseMatrix& mat, const MatrixZ& rhs)
      {
         rField.block(start, 0, mat.rows(), rField.cols()) = mat*rhs;
      }

      inline void setMatrixProduct(MatrixZ& rField, const int start, const SparseMatrixZ& mat, const Eigen::Ref<const MatrixZ>& rhs)
      {
         rField.block(start, 0, mat.rows(), rField.cols()) = mat*rhs;
      }

      inline void setMatrixProduct(DecoupledZMatrix& rField, const int start, const SparseMatrix& mat, const Eigen::Ref<const MatrixZ>& rhs)
      {
         int rows = mat.rows();
         int cols = rField.real().cols();
         rField.real().block(start, 0, rows, cols) = mat*rhs.real();
         rField.imag().block(start, 0, rows, cols) = mat*rhs.imag();
      }

      inline void setMatrixProduct(DecoupledZMatrix& rField, const int start, const SparseMatrix& mat, const Eigen::Ref<const Matrix>& rhsReal, const Eigen::Ref<const Matrix>& rhsImag)
      {
         assert(rField.real().rows() == rField.imag().rows());
         assert(rField.real().cols() == rField.imag().cols());

         int rows = mat.rows();
         int cols = rField.real().cols();
         rField.real().block(start, 0, rows, cols) = mat*rhsReal;
         rField.imag().block(start, 0, rows, cols) = mat*rhsImag;
      }

      inline void setMatrixProduct(DecoupledZMatrix& rField, const int start, const SparseMatrixZ& mat, const Eigen::Ref<const Matrix>& rhsReal, const Eigen::Ref<const Matrix>& rhsImag)
      {
         assert(rField.real().rows() == rField.imag().rows());
         assert(rField.real().cols() == rField.imag().cols());

         int rows = mat.rows();
         int cols = rField.real().cols();
         rField.real().block(start, 0, rows, cols) = mat.real()*rhsReal - mat.imag()*rhsImag;
         rField.imag().block(start, 0, rows, cols) = mat.real()*rhsImag + mat.imag()*rhsReal;
      }

      inline void setMatrixProduct(DecoupledZMatrix& rField, const int start, const SparseMatrixZ& mat, const Eigen::Ref<const MatrixZ>& rhs)
      {
         assert(rField.real().rows() == rField.imag().rows());
         assert(rField.real().cols() == rField.imag().cols());

         int rows = mat.rows();
         int cols = rField.real().cols();
         rField.real().block(start, 0, rows, cols) = mat.real()*rhs.real() - mat.imag()*rhs.imag();
         rField.imag().block(start, 0, rows, cols) = mat.real()*rhs.imag() + mat.imag()*rhs.real();
      }

      inline void addMatrixProduct(Matrix& rField, const int start, const SparseMatrixZ& mat, const Matrix& rhs)
      {
         throw Exception("Tried to use complex operator on real fields");
      }

      inline void addMatrixProduct(Matrix& rField, const int start, const SparseMatrix& mat, const Matrix& rhs)
      {
         rField.block(start, 0, mat.rows(), rField.cols()) += mat*rhs;
      }

      inline void addMatrixProduct(MatrixZ& rField, const int start, const SparseMatrix& mat, const MatrixZ& rhs)
      {
         rField.block(start, 0, mat.rows(), rField.cols()) += mat*rhs;
      }

      inline void addMatrixProduct(MatrixZ& rField, const int start, const SparseMatrixZ& mat, const Eigen::Ref<const MatrixZ>& rhs)
      {
         rField.block(start, 0, mat.rows(), rField.cols()) += mat*rhs;
      }

      inline void addMatrixProduct(DecoupledZMatrix& rField, const int start, const SparseMatrix& mat, const Eigen::Ref<const MatrixZ>& rhs)
      {
         int rows = mat.rows();
         int cols = rField.real().cols();
         rField.real().block(start, 0, rows, cols) += mat*rhs.real();
         rField.imag().block(start, 0, rows, cols) += mat*rhs.imag();
      }

      inline void addMatrixProduct(DecoupledZMatrix& rField, const int start, const SparseMatrix& mat, const Eigen::Ref<const Matrix>& rhsReal, const Eigen::Ref<const Matrix>& rhsImag)
      {
         assert(rField.real().rows() == rField.imag().rows());
         assert(rField.real().cols() == rField.imag().cols());

         int rows = mat.rows();
         int cols = rField.real().cols();
         rField.real().block(start, 0, rows, cols) += mat*rhsReal;
         rField.imag().block(start, 0, rows, cols) += mat*rhsImag;
      }

      inline void addMatrixProduct(DecoupledZMatrix& rField, const int start, const SparseMatrixZ& mat, const Eigen::Ref<const Matrix>& rhsReal, const Eigen::Ref<const Matrix>& rhsImag)
      {
         assert(rField.real().rows() == rField.imag().rows());
         assert(rField.real().cols() == rField.imag().cols());

         int rows = mat.rows();
         int cols = rField.real().cols();
         rField.real().block(start, 0, rows, cols) += mat.real()*rhsReal - mat.imag()*rhsImag;
         rField.imag().block(start, 0, rows, cols) += mat.real()*rhsImag + mat.imag()*rhsReal;
      }

      inline void addMatrixProduct(DecoupledZMatrix& rField, const int start, const SparseMatrixZ& mat, const Eigen::Ref<const MatrixZ>& rhs)
      {
         assert(rField.real().rows() == rField.imag().rows());
         assert(rField.real().cols() == rField.imag().cols());

         int rows = mat.rows();
         int cols = rField.real().cols();
         rField.real().block(start, 0, rows, cols) += mat.real()*rhs.real() - mat.imag()*rhs.imag();
         rField.imag().block(start, 0, rows, cols) += mat.real()*rhs.imag() + mat.imag()*rhs.real();
      }

      inline void setTopBlock(Matrix& rField, const int start, const int rows, const Matrix& rhs)
      {
         int cols = rField.cols();
         rField.block(start, 0, rows, cols) = rhs.topRows(rows);
      }

      inline void setTopBlock(MatrixZ& rField, const int start, const int rows, const MatrixZ& rhs)
      {
         int cols = rField.cols();
         rField.block(start, 0, rows, cols) = rhs.topRows(rows);
      }

      inline void setTopBlock(DecoupledZMatrix& rField, const int start, const int rows, const DecoupledZMatrix& rhs)
      {
         assert(rField.real().rows() == rField.imag().rows());
         assert(rField.real().cols() == rField.imag().cols());

         int cols = rField.real().cols();
         rField.real().block(start, 0, rows, cols) = rhs.real().topRows(rows);
         rField.imag().block(start, 0, rows, cols) = rhs.imag().topRows(rows);
      }
   }
}
}

#endif // MATRIXOPERATIONSINTERNAL_HPP
