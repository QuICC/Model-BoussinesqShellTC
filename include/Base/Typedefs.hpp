/** \file Typedefs.hpp
 *  \brief Definition of some useful typedefs used in the whole project
 *
 *  \mhdBug Needs test
 */

#ifndef TYPEDEFS_HPP
#define TYPEDEFS_HPP

/// Generate maximum precision diagnositics output
#define EIGEN_DEFAULT_IO_FORMAT IOFormat(Eigen::FullPrecision)

// Configuration includes
//
#include "SmartPointers/SharedPtrMacro.h"

// System includes
//
#include <complex>

// External includes
//
#include <Eigen/Core>
#include <Eigen/Sparse>

// Project includes
//

namespace GeoMHDiSCC {

   /**
    * @name Basic scalar types typedefs
    */
   //@{
   /// Typedef for floating point type value
   typedef double MHDFloat;
   /// Typedef for complex type value
   typedef std::complex<MHDFloat> MHDComplex;
   //@}

   /**
    * @name Fixed size array types typedefs
    */
   //@{
   /// Typedef for an array of boolean values
   typedef Eigen::Matrix<bool, 3, 1>   TriBool;
   //@}

   /**
    * @name Array types typedefs
    */
   //@{
   /// Typedef for an array of integer values
   typedef Eigen::Matrix<int, Eigen::Dynamic, 1>   ArrayI;
   /// Typedef for an array of float values
   typedef Eigen::Matrix<MHDFloat, Eigen::Dynamic, 1>   Array;
   /// Typedef for an array of complex values
   typedef Eigen::Matrix<MHDComplex, Eigen::Dynamic, 1>   ArrayZ;
   //@}

   /**
    * @name Matrix types typedefs
    */
   //@{
   /// Typedef for a matrix of float values
   typedef Eigen::Matrix<MHDFloat, Eigen::Dynamic, Eigen::Dynamic>   Matrix;
   /// Typedef for a matrix of complex values
   typedef Eigen::Matrix<MHDComplex, Eigen::Dynamic, Eigen::Dynamic>   MatrixZ;
   /// Typedef for a pair of real matrices used as a complex matrix
   typedef std::pair<Matrix,Matrix> DecoupledZMatrix;
   //@}

   /**
    * @name Sparse Matrix types typedefs
    */
   //@{
   /// Typedef for a sparse matrix of float values
   typedef Eigen::SparseMatrix<MHDFloat>   SparseMatrix;
   /// Typedef for a sparse matrix of complex values
   typedef Eigen::SparseMatrix<MHDComplex>   SparseMatrixZ;
   /// Typedef for a pair of real sparse matrices used as a sparse complex matrix
   typedef std::pair<SparseMatrix,SparseMatrix> DecoupledZSparse;
   //@}

   /**
    * @name Shared pointer typedefs
    */
   //@{
   /// Typedef for a smart reference counting pointer of an array of integers
   typedef SharedPtrMacro<ArrayI>   SharedArrayI;
   /// Typedef for an smart reference counting pointer of an array of reals
   typedef SharedPtrMacro<Array>   SharedArray;
   /// Typedef for an smart reference counting pointer of a matrix of reals
   typedef SharedPtrMacro<Matrix>   SharedMatrix;
   //@}
}

#endif // TYPEDEFS_HPP
