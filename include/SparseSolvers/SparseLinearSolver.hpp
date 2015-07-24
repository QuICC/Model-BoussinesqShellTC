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
#include "Base/MathConstants.hpp"
#include "Exceptions/Exception.hpp"
#include "Enums/SolveTiming.hpp"
#include "TypeSelectors/SparseSolverSelector.hpp"
#include "SparseSolvers/SparseSolverBase.hpp"
#include "SparseSolvers/SparseLinearSolverTools.hpp"

namespace GeoMHDiSCC {

namespace Solver {

   namespace internal {

      void addOperators(SparseMatrix& mat, const MHDFloat c, const DecoupledZSparse& decMat);

      void addOperators(SparseMatrixZ& mat, const MHDFloat c, const DecoupledZSparse& decMat);
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
          * @brief Initialise solution after data was copied
          */
         virtual void initSolutions();

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
          * @brief Prepare fields for implicit solve
          */
         virtual bool preSolve();

         /**
          * @brief Solve linear systems
          */
         void solve();

         /**
          * @brief Work on fields after implicit solve
          *
          * @param step    Current substep
          */
         virtual bool postSolve();

         /**
          * @brief Get the number of linear systems in solver
          */
         int nSystem() const;

         /**
          * @brief Set LHS matrix
          *
          * @param idx Index of the matrix
          */
         TOperator& rLHSMatrix(const MHDFloat id, const int idx);

         /**
          * @brief Add RHS and solution data storage
          * 
          * @param rows Number of rows of matrix
          * @param cols Number of columns required
          */
         virtual void addStorage(const int rows, const int cols);

         /**
          * @brief Build the scheme operators
          */
         void buildOperators(const int idx, const DecoupledZSparse& opA, const int size);

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

         /**
          * @brief Finished computation?
          */
         virtual bool finished();
         
      protected:
         /**
          * @brief Initialise the solver matrices storage
          *
          * @param n Size of matrices
          */
         void initMatrices(const MHDFloat id, const int n);

         /**
          * @brief Current ID of the solvers
          */
         MHDFloat mId;

         /**
          * @brief LHS operators
          */
         std::map<MHDFloat, std::vector<TOperator> >  mLHSMatrix;

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
         std::map<MHDFloat,std::vector<SharedPtrMacro<typename SparseSelector<TOperator>::Type > > >  mSolver;

      private:
   };

   template <typename TOperator,typename TData> SparseLinearSolver<TOperator,TData>::SparseLinearSolver(const int start, const SolveTiming::Id time)
      : SparseSolverBase(start, time), mId(0.0)
   {
   }

   template <typename TOperator,typename TData> SparseLinearSolver<TOperator,TData>::~SparseLinearSolver()
   {
   }

   template <typename TOperator,typename TData> bool SparseLinearSolver<TOperator,TData>::finished()
   {
      return true;
   }

   template <typename TOperator,typename TData> bool SparseLinearSolver<TOperator,TData>::preSolve()
   {
      return true;
   }

   template <typename TOperator,typename TData> void SparseLinearSolver<TOperator,TData>::solve()
   {
      // Set unused modes to zero
      for(int i = 0; i < this->mZeroIdx; ++i)
      {
         this->mSolution.at(i).setZero();
      }

      // Solve other modes
      typename std::map<MHDFloat, std::vector<SharedPtrMacro<typename SparseSelector<TOperator>::Type > > >::iterator sIt = this->mSolver.find(this->mId);
      for(size_t i = this->mZeroIdx; i < this->mRHSData.size(); i++)
      {
         #if defined GEOMHDISCC_MPI && defined GEOMHDISCC_MPISPSOLVE
            FrameworkMacro::syncSubComm(FrameworkMacro::SPECTRAL, i);
         #endif //define GEOMHDISCC_MPI && defined GEOMHDISCC_MPISPSOLVE

         //internal::solveWrapper<TOperator,TData>(this->mSolution.at(i), this->mSolver.at(i+start), this->mRHSData.at(i));
         internal::solveWrapper(this->mSolution.at(i), sIt->second.at(i), this->mRHSData.at(i));

         // Stop simulation if solve failed
         if(sIt->second.at(i)->info() != Eigen::Success)
         {
            throw Exception("Sparse direct solve failed!");
         }
      }
   }

   template <typename TOperator,typename TData> bool SparseLinearSolver<TOperator,TData>::postSolve()
   {
      return false;
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
      // Loop over matrices
      for(typename std::map<MHDFloat,std::vector<TOperator> >::iterator it = this->mLHSMatrix.begin(); it != this->mLHSMatrix.end(); ++it)
      {
         this->mSolver.insert(std::make_pair(it->first, std::vector<SharedPtrMacro<typename SparseSelector<TOperator>::Type > >()));

         typename std::map<MHDFloat, std::vector<SharedPtrMacro<typename SparseSelector<TOperator>::Type > > >::iterator sIt = this->mSolver.find(it->first);
         sIt->second.reserve(it->second.size());

         for(size_t i = 0; i < it->second.size(); ++i)
         {
            #if defined GEOMHDISCC_MPI && defined GEOMHDISCC_MPISPSOLVE
               FrameworkMacro::syncSubComm(FrameworkMacro::SPECTRAL, i);

               SharedPtrMacro<typename SparseSelector<TOperator>::Type >  solver(new typename SparseSelector<TOperator>::Type(FrameworkMacro::getSubComm(FrameworkMacro::SPECTRAL, i)));
            #else
               SharedPtrMacro<typename SparseSelector<TOperator>::Type >  solver(new typename SparseSelector<TOperator>::Type());
            #endif //define GEOMHDISCC_MPI && defined GEOMHDISCC_MPISPSOLVE

            sIt->second.push_back(solver);
         }
      }

      // Compute pattern and factorisation
      this->updateSolver();
   }

   template <typename TOperator,typename TData> void SparseLinearSolver<TOperator,TData>::updateSolver()
   {
      for(typename std::map<MHDFloat,std::vector<TOperator> >::iterator it = this->mLHSMatrix.begin(); it != this->mLHSMatrix.end(); ++it)
      {
         // Compute factorisation
         for(size_t i = 0; i < it->second.size(); i++)
         {
            if(static_cast<int>(i) % this->nSystem() >= this->mZeroIdx)
            {
               #if defined GEOMHDISCC_MPI && defined GEOMHDISCC_MPISPSOLVE
                  FrameworkMacro::syncSubComm(FrameworkMacro::SPECTRAL, i);
               #endif //define GEOMHDISCC_MPI && defined GEOMHDISCC_MPISPSOLVE

               typename std::map<MHDFloat, std::vector<SharedPtrMacro<typename SparseSelector<TOperator>::Type > > >::iterator sIt = this->mSolver.find(it->first);
               // Safety assert to make sur matrix is compressed
               assert(it->second.at(i).isCompressed());

               sIt->second.at(i)->compute(it->second.at(i));

               // Stop simulation if factorization failed
               if(sIt->second.at(i)->info() != Eigen::Success)
               {
                  throw Exception("Matrix factorization failed!");
               }
            }
         }
      }
   }

   template <typename TOperator,typename TData> void SparseLinearSolver<TOperator,TData>::initMatrices(const int n)
   {
      this->initMatrices(0, n);
   }

   template <typename TOperator,typename TData> void SparseLinearSolver<TOperator,TData>::initMatrices(const MHDFloat id, const int n)
   {
      // Do not reinitialise if work already done by other field
      if(this->mLHSMatrix.count(id) == 0 || this->mLHSMatrix.find(id)->second.size() == 0)
      {  
         if(this->mLHSMatrix.count(id) == 0)
         {
            this->mLHSMatrix.insert(std::make_pair(id, std::vector<TOperator>()));
         }

         typename std::map<MHDFloat,std::vector<TOperator> >::iterator it = this->mLHSMatrix.find(id);

         // Reserve space for the LHS matrices
         it->second.reserve(n);

         // Initialise storage
         for(int i = 0; i < n; ++i)
         {
            // Create storage for LHS matrices
            it->second.push_back(TOperator());
         }
      }
   }

   template <typename TOperator,typename TData> void SparseLinearSolver<TOperator,TData>::initSolutions()
   {
      // Nothing to be done in general.
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

   template <typename TOperator,typename TData> void SparseLinearSolver<TOperator,TData>::buildOperators(const int idx, const DecoupledZSparse& opA, const int size)
   {
      this->rLHSMatrix(0, idx).resize(size, size);
      Solver::internal::addOperators(this->rLHSMatrix(0, idx), 1.0, opA);
   }

   template <typename TOperator,typename TData> int SparseLinearSolver<TOperator,TData>::nSystem() const
   {
      return this->mRHSData.size();
   }

   template <typename TOperator,typename TData> TOperator& SparseLinearSolver<TOperator,TData>::rLHSMatrix(const MHDFloat id, const int idx)
   {
      return this->mLHSMatrix.find(id)->second.at(idx);
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

   namespace internal {

      inline void addOperators(SparseMatrix& mat, const MHDFloat c, const DecoupledZSparse& decMat)
      {
         assert(decMat.real().rows() > 0);
         assert(decMat.real().cols() > 0);
         assert(decMat.imag().size() == 0 || decMat.imag().nonZeros() == 0);

         if(c != 1.0)
         {
            mat += c*decMat.real();
         } else
         {
            mat += decMat.real();
         }
      }

      inline void addOperators(SparseMatrixZ& mat, const MHDFloat c, const DecoupledZSparse& decMat)
      {
         assert(decMat.real().rows() > 0);
         assert(decMat.real().cols() > 0);
         assert(decMat.imag().rows() > 0);
         assert(decMat.imag().cols() > 0);
         assert(decMat.real().rows() == decMat.imag().rows());
         assert(decMat.real().cols() == decMat.imag().cols());

         if(c != 1.0)
         {
            mat += c*decMat.real().cast<MHDComplex>() + c*Math::cI*decMat.imag();
         } else
         {
            mat += decMat.real().cast<MHDComplex>() + Math::cI*decMat.imag();
         }
      }
   }
}
}

#endif // SPARSELINEARSOLVER_HPP
