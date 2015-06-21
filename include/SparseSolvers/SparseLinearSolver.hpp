/** 
 * @file SparseLinearSolver.hpp
 * @brief Implementation of a templated (coupled) linear solver structure
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef SPARSELINEARSOLVER_HPP
#define SPARSELINEARSOLVER_HPP

// Configuration includes
//
#include "Framework/FrameworkMacro.h"
#include "SmartPointers/SharedPtrMacro.h"

// System includes
//

// External includes
//

// Project includes
//
#include "Enums/SolveTiming.hpp"
#include "TypeSelectors/SparseSolverSelector.hpp"
#include "SparseSolvers/SparseSolverBase.hpp"
#include "SparseSolvers/SparseLinearSolverTools.hpp"

namespace GeoMHDiSCC {

namespace Solver {

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
          * @param time    Solver timing with respect to timestepping
          */
         SparseLinearSolver(const int start, const SolveTiming::Id time);

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
          * @brief Set solver RHS data to zero
          */
         void zeroSolver();

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
         std::vector<SharedPtrMacro<typename SparseSelector<TOperator>::Type > >  mSolver;

      private:
   };

   template <typename TOperator,typename TData> SparseLinearSolver<TOperator,TData>::SparseLinearSolver(const int start, const SolveTiming::Id time)
      : SparseSolverBase(start, time)
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
         #if defined GEOMHDISCC_MPI && defined GEOMHDISCC_MPISPSOLVE
            FrameworkMacro::syncSubComm(FrameworkMacro::SPECTRAL, i);
         #endif //define GEOMHDISCC_MPI && defined GEOMHDISCC_MPISPSOLVE

         //internal::solveWrapper<TOperator,TData>(this->mSolution.at(i), this->mSolver.at(i+start), this->mRHSData.at(i));
         internal::solveWrapper(this->mSolution.at(i), this->mSolver.at(i+start), this->mRHSData.at(i));

         // Stop simulation if solve failed
         if(this->mSolver.at(i+start)->info() != Eigen::Success)
         {
            throw Exception("Sparse direct solve failed!");
         }
      }
   }

   template <typename TOperator,typename TData> void SparseLinearSolver<TOperator,TData>::zeroSolver()
   {
      // Set solver RHS to zero
      for(int i = 0; i < this->mRHSData.size(); ++i)
      {
         this->mRHSData.at(i).setZero();
      }
   }

   template <typename TOperator,typename TData> void SparseLinearSolver<TOperator,TData>::initSolver()
   {
      // Initialise solver
      this->mSolver.reserve(this->mLHSMatrix.size());

      for(size_t i = 0; i < this->mLHSMatrix.size(); i++)
      {
         #if defined GEOMHDISCC_MPI && defined GEOMHDISCC_MPISPSOLVE
            FrameworkMacro::syncSubComm(FrameworkMacro::SPECTRAL, i % this->nSystem());

            SharedPtrMacro<typename SparseSelector<TOperator>::Type >  solver(new typename SparseSelector<TOperator>::Type(FrameworkMacro::getSubComm(FrameworkMacro::SPECTRAL, i % this->nSystem())));
         #else
            SharedPtrMacro<typename SparseSelector<TOperator>::Type >  solver(new typename SparseSelector<TOperator>::Type());
         #endif //define GEOMHDISCC_MPI && defined GEOMHDISCC_MPISPSOLVE

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
            #if defined GEOMHDISCC_MPI && defined GEOMHDISCC_MPISPSOLVE
               FrameworkMacro::syncSubComm(FrameworkMacro::SPECTRAL, i % this->nSystem());
            #endif //define GEOMHDISCC_MPI && defined GEOMHDISCC_MPISPSOLVE

            // Safety assert to make sur matrix is compressed
            assert(this->mLHSMatrix.at(i).isCompressed());

            this->mSolver.at(i)->compute(this->mLHSMatrix.at(i));

            // Stop simulation if factorization failed
            if(this->mSolver.at(i)->info() != Eigen::Success)
            {
               throw Exception("Matrix factorization failed!");
            }
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
}
}

#endif // SPARSELINEARSOLVER_HPP
