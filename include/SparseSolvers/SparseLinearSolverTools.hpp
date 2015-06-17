/** 
 * @file SparseLinearSolverTools.hpp
 * @brief Implementation of a couple of templated linear solver helper functions
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef SPARSELINEARSOLVERTOOLS_HPP
#define SPARSELINEARSOLVERTOOLS_HPP

// Configuration includes
//
#include "SmartPointers/SharedPtrMacro.h"

// System includes
//

// External includes
//

// Project includes
//
#include "TypeSelectors/SparseSolverSelector.hpp"
#include "Base/Typedefs.hpp"

namespace GeoMHDiSCC {

namespace Solver {

namespace internal {

   /**
    * @brief Simpler wrapper around solver solve call
    */
   void solveWrapper(Matrix& rSolution, SparseSelector<SparseMatrix>::Type& solver, const Matrix& rhs);
   void solveWrapper(MatrixZ& rSolution, SparseSelector<SparseMatrix>::Type& solver, const MatrixZ& rhs);
//   void solveWrapper(MatrixZ& rSolution, SparseSelector<SparseMatrixZ>::Type& solver, const MatrixZ& rhs);
   void solveWrapper(DecoupledZMatrix& rSolution, SparseSelector<SparseMatrix>::Type& solver, const DecoupledZMatrix& rhs);

   /**
    * @brief Simpler wrapper around solver solve call for shared pointer to solver
    */
   void solveWrapper(Matrix& rSolution, SharedPtrMacro< SparseSelector<SparseMatrix>::Type > solver, const Matrix& rhs);
//   void solveWrapper(MatrixZ& rSolution, SharedPtrMacro< SparseSelector<SparseMatrix>::Type > solver, const MatrixZ& rhs);
   void solveWrapper(MatrixZ& rSolution, SharedPtrMacro< SparseSelector<SparseMatrixZ>::Type > solver, const MatrixZ& rhs);
   void solveWrapper(DecoupledZMatrix& rSolution, SharedPtrMacro< SparseSelector<SparseMatrix>::Type > solver, const DecoupledZMatrix& rhs);
}

namespace internal {

   inline void solveWrapper(Matrix& rSolution, SparseSelector<SparseMatrix>::Type& solver, const Matrix& rhs)
   {
      rSolution = solver.solve(rhs);

      // Safety assert for successful solve
      assert(solver.info() == Eigen::Success);
   }

   inline void solveWrapper(MatrixZ& rSolution, SparseSelector<SparseMatrix>::Type& solver, const MatrixZ& rhs)
   {
      Matrix tmpIn(rhs.rows(), rhs.cols());
      Matrix tmpOut(rhs.rows(), rhs.cols());

      tmpIn = rhs.real();
      tmpOut = solver.solve(tmpIn);
      rSolution.real() = tmpOut;

      // Safety assert for successful solve
      assert(solver.info() == Eigen::Success);

      tmpIn = rhs.imag();
      tmpOut = solver.solve(tmpIn);
      rSolution.imag() = tmpOut;

      // Safety assert for successful solve
      assert(solver.info() == Eigen::Success);
   }

//   inline void solveWrapper(MatrixZ& rSolution, SparseSelector<SparseMatrixZ>::Type& solver, const MatrixZ& rhs)
//   {
//      rSolution = solver.solve(rhs);
//
//      // Safety assert for successful solve
//      assert(solver.info() == Eigen::Success);
//   }

   inline void solveWrapper(DecoupledZMatrix& rSolution, SparseSelector<SparseMatrix>::Type& solver, const DecoupledZMatrix& rhs)
   {
      rSolution.real() = solver.solve(rhs.real());

      // Safety assert for successful solve
      assert(solver.info() == Eigen::Success);

      rSolution.imag() = solver.solve(rhs.imag());

      // Safety assert for successful solve
      assert(solver.info() == Eigen::Success);
   }

   inline void solveWrapper(Matrix& rSolution, SharedPtrMacro< SparseSelector<SparseMatrix>::Type > spSolver, const Matrix& rhs)
   {
      rSolution = spSolver->solve(rhs);

      // Safety assert for successful solve
      assert(spSolver->info() == Eigen::Success);
   }

//   inline void solveWrapper(MatrixZ& rSolution, SharedPtrMacro< SparseSelector<SparseMatrix>::Type > spSolver, const MatrixZ& rhs)
//   {
//      Matrix tmp(rhs.rows(), rhs.cols());
//      tmp = rhs.real();
//      rSolution.real() = spSolver->solve(tmp);
//
//      // Safety assert for successful solve
//      assert(spSolver->info() == Eigen::Success);
//
//      tmp = rhs.imag();
//      rSolution.imag() = spSolver->solve(tmp);
//
//      // Safety assert for successful solve
//      assert(spSolver->info() == Eigen::Success);
//   }

   inline void solveWrapper(MatrixZ& rSolution, SharedPtrMacro< SparseSelector<SparseMatrixZ>::Type > spSolver, const MatrixZ& rhs)
   {
      rSolution = spSolver->solve(rhs);

      // Safety assert for successful solve
      assert(spSolver->info() == Eigen::Success);
   }

   inline void solveWrapper(DecoupledZMatrix& rSolution, SharedPtrMacro<SparseSelector<SparseMatrix>::Type > spSolver, const DecoupledZMatrix& rhs)
   {
      rSolution.real() = spSolver->solve(rhs.real());

      // Safety assert for successful solve
      assert(spSolver->info() == Eigen::Success);

      rSolution.imag() = spSolver->solve(rhs.imag());

      // Safety assert for successful solve
      assert(spSolver->info() == Eigen::Success);
   }
}
}
}

#endif // SPARSELINEARSOLVERTOOLS_HPP
