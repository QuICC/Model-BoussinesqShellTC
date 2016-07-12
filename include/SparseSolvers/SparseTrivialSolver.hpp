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
#include "Enums/SolveTiming.hpp"
#include "SparseSolvers/SparseSolverBase.hpp"

namespace GeoMHDiSCC {

namespace Solver {

   /**
    * @brief Implementation of a templated (coupled) trivial solver structure
    */
   template <typename TOperator,typename TData> class SparseTrivialSolver: public SparseSolverBase
   {
      public:
         /**
          * @brief Constructor
          *
          * @param start   Starting index (for example without m=0)
          * @param time    Solver timing with respect to timestepping
          * @param expTime Explicit linear timing with respect to nonlinear calculation
          */
         SparseTrivialSolver(const int start, const SolveTiming::Id time);

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
          * @brief Initialise solution after data was copied
          */
         virtual void initSolutions();

         /**
          * @brief Set solver RHS data to zero
          */
         void zeroSolver();

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

   template <typename TOperator,typename TData> SparseTrivialSolver<TOperator,TData>::SparseTrivialSolver(const int start, const SolveTiming::Id time)
      : SparseSolverBase(start, time)
   {
      this->setInitialized();
   }

   template <typename TOperator,typename TData> SparseTrivialSolver<TOperator,TData>::~SparseTrivialSolver()
   {
   }

   template <typename TOperator,typename TData> void SparseTrivialSolver<TOperator,TData>::addStorage(const int rows, const int cols)
   {
      // Assert for non zero rows and columns
      assert(rows > 0);
      assert(cols > 0);

      // Add storage for solution
      this->mSolution.push_back(TData(rows,cols));
      this->mSolution.back().setZero();
   }

   template <typename TOperator,typename TData> void SparseTrivialSolver<TOperator,TData>::initSolutions()
   {
   }

   template <typename TOperator,typename TData> void SparseTrivialSolver<TOperator,TData>::zeroSolver()
   {
      // Set solver RHS to zero
      for(unsigned int i = 0; i < this->mSolution.size(); ++i)
      {
         this->mSolution.at(i).setZero();
      }
   }

   template <typename TOperator,typename TData> int SparseTrivialSolver<TOperator,TData>::nSystem() const
   {
      return this->mSolution.size();
   }

   template <typename TOperator,typename TData> const TData& SparseTrivialSolver<TOperator,TData>::solution(const int idx) const
   {
      return this->mSolution.at(idx);
   }

   template <typename TOperator,typename TData> TData& SparseTrivialSolver<TOperator,TData>::rSolution(const int idx)
   {
      return this->mSolution.at(idx);
   }

   template <typename TOperator,typename TData> TData& SparseTrivialSolver<TOperator,TData>::rRHSData(const int idx)
   {
      // WARNING: this is the same as rSolution. It's used to simplify some implementations
      return this->mSolution.at(idx);
   }
}
}

#endif // SPARSETRIVIALSOLVER_HPP
