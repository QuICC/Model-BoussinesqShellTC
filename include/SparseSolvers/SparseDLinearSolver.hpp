/** \file SparseDLinearSolver.hpp
 *  \brief Implementation of a real (coupled) linear solver structure
 */

#ifndef SPARSEDLINEARSOLVER_HPP
#define SPARSEDLINEARSOLVER_HPP

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
#include "SparseSolvers/SparseLinearSolverBase.hpp"

namespace GeoMHDiSCC {

namespace Solver {

   /**
    *  \brief Implementation of a real (coupled) linear solver structure
    */
   class SparseDLinearSolver: public SparseLinearSolverBase
   {
      public:
         /**
          * @brief Constructor
          *
          * @param start   Starting index (for example without m=0)
          */
         SparseDLinearSolver(const int start);

         /**
          * @brief Destructor
          */
         virtual ~SparseDLinearSolver();

         /**
          * @brief Initialise the solver matrices storage
          *
          * @param n Size of matrices
          */
         void initMatrices(const int n);

         /**
          * @brief Initialise solver
          */
         void initSolver();

         /**
          * @brief Update solver
          */
         void updateSolver();

         /**
          * @brief Compute the RHS of the linear systems
          *
          * @param step    Current substep
          */
         void computeRHS(const int step);

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
         SparseMatrix& rLHSMatrix(const int idx);

         /**
          * @brief Add RHS and solution data storage
          * 
          * @param rows Number of rows of matrix
          * @param cols Number of columns required
          */
         void addStorage(const int rows, const int cols);

         /**
          * @brief Set RHS data
          *
          * @param idx   Index of the data
          */
         DecoupledZMatrix& rRHSData(const int idx);

         /**
          * @brief Get solution data
          *
          * @param idx   Index of the data
          */
         const DecoupledZMatrix& solution(const int idx) const;

         /**
          * @brief Set solution data
          *
          * @param idx   Index of the data
          */
         DecoupledZMatrix& rSolution(const int idx);
         
      protected:
         /**
          * @brief Real LHS operators of the timestepped equations
          */
         std::vector<SparseMatrix>   mLHSMatrix;

         /**
          * @brief Storage for linear solve's RHS
          */
         std::vector<DecoupledZMatrix>  mRHSData;

         /**
          * @brief Storage for old nonlinear RHS
          */
         std::vector<DecoupledZMatrix>  mRHSOld;

         /**
          * @brief Storage for solution of linear solve
          */
         std::vector<DecoupledZMatrix>  mSolution;

         /**
          * @brief Create sparse solvers
          */
         std::vector<SharedPtrMacro<SparseSolverMacro<SparseMatrix> > >  mSolver;

      private:
   };

   inline int SparseDLinearSolver::nSystem() const
   {
      return this->mRHSData.size();
   }

   inline DecoupledZMatrix& SparseDLinearSolver::rRHSData(const int idx)
   {
      return this->mRHSData.at(idx);
   }

   inline const DecoupledZMatrix& SparseDLinearSolver::solution(const int idx) const
   {
      return this->mSolution.at(idx);
   }

   inline DecoupledZMatrix& SparseDLinearSolver::rSolution(const int idx)
   {
      return this->mSolution.at(idx);
   }

   /// Typedef for a shared pointer of a SparseDLinearSolver
   typedef SharedPtrMacro<SparseDLinearSolver>  SharedSparseDLinearSolver;
}
}

#endif // SPARSEDLINEARSOLVER_HPP
