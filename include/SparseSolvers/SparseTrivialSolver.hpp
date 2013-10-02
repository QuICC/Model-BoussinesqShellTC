/** 
 * @file SparseTrivialSolver.hpp
 * @brief Implementation of a templated (coupled) trivial solver structure
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef SPARSETRIVIALSOLVER_HPP
#define SPARSETRIVIALSOLVER_HPP

// Configuration includes
//
#include "SmartPointers/SharedPtrMacro.h"

// System includes
//

// External includes
//

// Project includes
//
#include "SparseSolvers/SparseSolverBase.hpp"

namespace GeoMHDiSCC {

namespace Solver {

   /**
    * @brief Implementation of a templated (coupled) trivial solver structure
    */
   template <typename TData> class SparseTrivialSolver: public SparseSolverBase
   {
      public:
         /**
          * @brief Constructor
          *
          * @param start   Starting index (for example without m=0)
          */
         SparseTrivialSolver(const int start);

         /**
          * @brief Destructor
          */
         virtual ~SparseTrivialSolver();

         /**
          * @brief Get the number of systems in solver
          */
         int nSystem() const;

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
          * @brief Storage for solution of linear solve
          */
         std::vector<TData>  mSolution;

      private:
   };

   /// Typedef for a real data real matrix trivial solver
   typedef SparseTrivialSolver<Matrix> SparseRTrivialSolver;
   /// Typedef for a shared real data real matrix trivial solver
   typedef SharedPtrMacro<SparseRTrivialSolver> SharedSparseRTrivialSolver;
   /// Typedef for a complex data real matrix trivial solver
   typedef SparseTrivialSolver<DecoupledZMatrix> SparseRZTrivialSolver;
   /// Typedef for a shared real data real matrix trivial solver
   typedef SharedPtrMacro<SparseRZTrivialSolver> SharedSparseRZTrivialSolver;
   /// Typedef for a complex data complex matrix trivial solver
   typedef SparseTrivialSolver<MatrixZ> SparseZTrivialSolver;
   /// Typedef for a shared real data real matrix trivial solver
   typedef SharedPtrMacro<SparseZTrivialSolver> SharedSparseZTrivialSolver;

   template <typename TData> SparseTrivialSolver<TData>::SparseTrivialSolver(const int start)
      : SparseSolverBase(start)
   {
   }

   template <typename TData> SparseTrivialSolver<TData>::~SparseTrivialSolver()
   {
   }

   template <typename TData> void SparseTrivialSolver<TData>::addStorage(const int rows, const int cols)
   {
      // Assert for non zero rows and columns
      assert(rows > 0);
      assert(cols > 0);

      // Add storage for solution
      this->mSolution.push_back(TData(rows,cols));
      this->mSolution.back().setZero();
   }

   template <typename TData> int SparseTrivialSolver<TData>::nSystem() const
   {
      return this->mSolution.size();
   }

   template <typename TData> const TData& SparseTrivialSolver<TData>::solution(const int idx) const
   {
      return this->mSolution.at(idx);
   }

   template <typename TData> TData& SparseTrivialSolver<TData>::rSolution(const int idx)
   {
      return this->mSolution.at(idx);
   }

   template <typename TData> TData& SparseTrivialSolver<TData>::rRHSData(const int idx)
   {
      // WARNING: this is the same as rSolution. It's used to simplify some implementations
      return this->mSolution.at(idx);
   }
}
}

#endif // SPARSETRIVIALSOLVER_HPP
