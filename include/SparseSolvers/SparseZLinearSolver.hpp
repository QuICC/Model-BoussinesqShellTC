/** \file SparseZLinearSolver.hpp
 *  \brief Implementation of a complex (coupled) linear solver structure
 */

#ifndef SPARSEZLINEARSOLVER_HPP
#define SPARSEZLINEARSOLVER_HPP

// Configuration includes
//
#include "SmartPointers/SharedPtrMacro.h"

// System includes
//

// External includes
//

// Project includes
//
#include "SparseSolvers/SparseSolverMacro.h"
#include "SparseSolvers/SparseSolverBase.hpp"

namespace GeoMHDiSCC {

namespace Solver {

   /**
    *  \brief Implementation of a complex (coupled) linear solver structure
    */
   class SparseZLinearSolver: public SparseSolverBase
   {
      public:
         /**
          * @brief Constructor
          *
          * @param start   Starting index (for example without m=0)
          */
         SparseZLinearSolver(const int start);

         /**
          * @brief Destructor
          */
         virtual ~SparseZLinearSolver();

         /**
          * @brief Initialise the solver matrices storage
          *
          * @param n Size of matrices
          */
         virtual void initMatrices(const int n);

         /**
          * @brief Initialise solver
          */
         void initSolver();

         /**
          * @brief Update solver
          */
         void updateSolver();

         /**
          * @brief Solve linear systems
          *
          * @param step    Current substep
          */
         void solve(const int step);

         /**
          * @brief Get the number of linear systems in solver
          */
         int nSystem() const;

         /**
          * @brief Set LHS matrix
          *
          * @param idx Index of the matrix
          */
         SparseMatrixZ& rLHSMatrix(const int idx);

         /**
          * @brief Add RHS and solution data storage
          * 
          * @param rows Number of rows of matrix
          * @param cols Number of columns required
          */
         virtual void addStorage(const int rows, const int cols);

         /**
          * @brief Set RHS data
          *
          * @param idx   Index of the data
          */
         MatrixZ& rRHSData(const int idx);

         /**
          * @brief Get solution data
          *
          * @param idx   Index of the data
          */
         const MatrixZ& solution(const int idx) const;

         /**
          * @brief Set solution data
          *
          * @param idx   Index of the data
          */
         MatrixZ& rSolution(const int idx);
         
      protected:
         /**
          * @brief Complex LHS operators of the timestepped equations
          */
         std::vector<SparseMatrixZ>   mLHSMatrix;

         /**
          * @brief Storage for linear solve's RHS
          */
         std::vector<MatrixZ>  mRHSData;

         /**
          * @brief Storage for solution of linear solve
          */
         std::vector<MatrixZ>  mSolution;

         /**
          * @brief Create sparse solvers
          */
         std::vector<SharedPtrMacro<SparseSolverMacro<SparseMatrixZ> > >  mSolver;

      private:
   };

   /// Typedef for a shared pointer of a SparseZLinearSolver
   typedef SharedPtrMacro<SparseZLinearSolver>  SharedSparseZLinearSolver;
}
}

#endif // SPARSEZLINEARSOLVER_HPP
