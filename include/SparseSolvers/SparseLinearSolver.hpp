/** 
 * @file SparseLinearSolver.hpp
 * @brief Implementation of a templated (coupled) linear solver structure
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef SPARSELINEARSOLVER_HPP
#define SPARSELINEARSOLVER_HPP

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
#include "SparseSolvers/SparseSolverBase.hpp"

namespace GeoMHDiSCC {

namespace Solver {

   namespace internal
   {
      template <typename TOperator,typename TData> inline void solveWrapper(TData& rSolution, SharedPtrMacro<typename SparseSolverSelector<TOperator>::SolverType > solver, const TData& rhs);
   }

   /**
    * @brief Implementation of a templated (coupled) linear solver structure
    */
   template <typename TOperator, typename TData> class SparseLinearSolver: public SparseSolverBase
   {
      public:
         /**
          * @brief Constructor
          *
          * @param start   Starting index (for example without m=0)
          */
         SparseLinearSolver(const int start);

         /**
          * @brief Destructor
          */
         virtual ~SparseLinearSolver();

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
         TOperator& rLHSMatrix(const int idx);

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
         TData& rRHSData(const int idx);

         /**
          * @brief Get solution data
          *
          * @param idx   Index of the data
          */
         const TData& solution(const int idx) const;

         /**
          * @brief Set solution data
          *
          * @param idx   Index of the data
          */
         TData& rSolution(const int idx);
         
      protected:
         /**
          * @brief Complex LHS operators of the timestepped equations
          */
         std::vector<TOperator>   mLHSMatrix;

         /**
          * @brief Storage for linear solve's RHS
          */
         std::vector<TData>  mRHSData;

         /**
          * @brief Storage for solution of linear solve
          */
         std::vector<TData>  mSolution;

         /**
          * @brief Create sparse solvers
          */
         std::vector<SharedPtrMacro<typename SparseSolverSelector<TOperator>::SolverType > >  mSolver;

      private:
   };

   /// Typedef for a real data real matrix solver
   typedef SparseLinearSolver<SparseMatrix,Matrix> SparseRLinearSolver;
   /// Typedef for a shared real data real matrix solver
   typedef SharedPtrMacro<SparseRLinearSolver> SharedSparseRLinearSolver;
   /// Typedef for a complex data real matrix solver
   typedef SparseLinearSolver<SparseMatrix,DecoupledZMatrix> SparseRZLinearSolver;
   /// Typedef for a shared real data real matrix solver
   typedef SharedPtrMacro<SparseRZLinearSolver> SharedSparseRZLinearSolver;
   /// Typedef for a complex data complex matrix solver
   typedef SparseLinearSolver<SparseMatrixZ,MatrixZ> SparseZLinearSolver;
   /// Typedef for a shared real data real matrix solver
   typedef SharedPtrMacro<SparseZLinearSolver> SharedSparseZLinearSolver;

   template <typename TOperator,typename TData> SparseLinearSolver<TOperator,TData>::SparseLinearSolver(const int start)
      : SparseSolverBase(start)
   {
   }

   template <typename TOperator,typename TData> SparseLinearSolver<TOperator,TData>::~SparseLinearSolver()
   {
   }

   template <typename TOperator,typename TData> void SparseLinearSolver<TOperator,TData>::solve(const int step)
   {
      int start = step*this->nSystem();

      // Set unused modes to zero
      for(int i = 0; i < this->mZeroIdx; ++i)
      {
         this->mSolution.at(i).setZero();
      }

      // Solve other modes
      for(size_t i = this->mZeroIdx; i < this->mRHSData.size(); i++)
      {
         internal::solveWrapper<TOperator,TData>(this->mSolution.at(i), this->mSolver.at(i+start), this->mRHSData.at(i));
      }
   }

   template <typename TOperator,typename TData> void SparseLinearSolver<TOperator,TData>::initSolver()
   {
      // Initialise solver
      this->mSolver.reserve(this->mLHSMatrix.size());
      for(size_t i = 0; i < this->mLHSMatrix.size(); i++)
      {
         SharedPtrMacro<typename SparseSolverSelector<TOperator>::SolverType >  solver(new typename SparseSolverSelector<TOperator>::SolverType());

         this->mSolver.push_back(solver);
      }

      // Compute pattern and factorisation
      this->updateSolver();
   }

   template <typename TOperator,typename TData> void SparseLinearSolver<TOperator,TData>::updateSolver()
   {
      // Compute factorisation
      for(size_t i = 0; i < this->mLHSMatrix.size(); i++)
      {
         if(static_cast<int>(i) % this->nSystem() >= this->mZeroIdx)
         {
            // Safety assert to make sur matrix is compressed
            assert(this->mLHSMatrix.at(i).isCompressed());

            this->mSolver.at(i)->compute(this->mLHSMatrix.at(i));

            // Safety assert for successful factorisation
            assert(this->mSolver.at(i)->info() == Eigen::Success);
         }
      }
   }

   template <typename TOperator,typename TData> void SparseLinearSolver<TOperator,TData>::initMatrices(const int n)
   {
      // Do not reinitialise if work already done by other field
      if(this->mLHSMatrix.size() == 0)
      {
         // Reserve space for the LHS matrices
         this->mLHSMatrix.reserve(n);

         // Initialise storage
         for(int i = 0; i < n; ++i)
         {
            // Create storage for LHS matrices
            this->mLHSMatrix.push_back(TOperator());
         }
      }
   }

   template <typename TOperator,typename TData> void SparseLinearSolver<TOperator,TData>::addStorage(const int rows, const int cols)
   {
      // Assert for non zero rows and columns
      assert(rows > 0);
      assert(cols > 0);

      // Add storage for RHS data
      this->mRHSData.push_back(TData(rows,cols));
      this->mRHSData.back().setZero();

      // Add storage for solution
      this->mSolution.push_back(TData(rows,cols));
      this->mSolution.back().setZero();
   }

   template <typename TOperator,typename TData> int SparseLinearSolver<TOperator,TData>::nSystem() const
   {
      return this->mRHSData.size();
   }

   template <typename TOperator,typename TData> TOperator& SparseLinearSolver<TOperator,TData>::rLHSMatrix(const int idx)
   {
      return this->mLHSMatrix.at(idx);
   }

   template <typename TOperator,typename TData> TData& SparseLinearSolver<TOperator,TData>::rRHSData(const int idx)
   {
      return this->mRHSData.at(idx);
   }

   template <typename TOperator,typename TData> const TData& SparseLinearSolver<TOperator,TData>::solution(const int idx) const
   {
      return this->mSolution.at(idx);
   }

   template <typename TOperator,typename TData> TData& SparseLinearSolver<TOperator,TData>::rSolution(const int idx)
   {
      return this->mSolution.at(idx);
   }

   namespace internal
   {
   template <typename TOperator,typename TData> inline void solveWrapper(TData& rSolution, SharedPtrMacro<typename SparseSolverSelector<TOperator>::SolverType > solver, const TData& rhs)
   {
      rSolution = solver->solve(rhs);

      // Safety assert for successful solve
      assert(solver->info() == Eigen::Success);
   }

   template <> inline void solveWrapper<SparseMatrix,DecoupledZMatrix>(DecoupledZMatrix& rSolution, SharedPtrMacro<SparseSolverSelector<SparseMatrix>::SolverType > solver, const DecoupledZMatrix& rhs)
   {
      rSolution.real() = solver->solve(rhs.real());

      // Safety assert for successful solve
      assert(solver->info() == Eigen::Success);

      rSolution.imag() = solver->solve(rhs.imag());

      // Safety assert for successful solve
      assert(solver->info() == Eigen::Success);
   }
   }
}
}

#endif // SPARSELINEARSOLVER_HPP