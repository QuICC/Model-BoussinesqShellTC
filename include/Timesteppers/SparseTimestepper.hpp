/** 
 * @file SparseTimestepper.hpp
 * @brief Implementation of a templated (coupled) equation timestepper
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef SPARSETIMESTEPPER_HPP
#define SPARSETIMESTEPPER_HPP

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
#include "TypeSelectors/TimeSchemeSelector.hpp"
#include "SparseSolvers/SparseLinearSolver.hpp"

namespace GeoMHDiSCC {

namespace Timestep {

   namespace internal
   {
      template <typename TOperator,typename TData> void computeRHSNoMemory(const int step, TData& rRHS, const TOperator& mat, const TData& sol, TData& rOld);

      template <typename TOperator,typename TData> void computeRHSMemory(const int step, TData& rRHS, const TOperator& mat, const TData& sol, TData& rOld, TData& rTmp);
   }

   /**
    * @brief Implementation of a templated (coupled) equation timestepper
    */
   template <typename TOperator,typename TData> class SparseTimestepper: public Solver::SparseLinearSolver<TOperator,TData>
   {
      public:
         /**
          * @brief Constructor
          *
          * @param start   Starting index (for example without m=0)
          */
         SparseTimestepper(const int start);

         /**
          * @brief Destructor
          */
         virtual ~SparseTimestepper();

         /**
          * @brief Initialise the solver matrices storage
          *
          * @param n Size of matrices
          */
         virtual void initMatrices(const int n);

         /**
          * @brief Compute the RHS of the linear systems
          *
          * @param step    Current substep
          */
         void computeRHS(const int step);

         /**
          * @brief Update the LHS matrix with new timedependence
          *
          * @param lhsCoeff   New coefficient for LHS time dependent part
          * @param rhsCoeff   New coefficient for RHS time dependent part
          * @param step       Timestep scheme substep
          */
         void updateTimeMatrix(const MHDFloat lhsCoeff, const MHDFloat rhsCoeff, const int step);

         /**
          * @brief Set RHS matrix
          *
          * @param idx Index of the matrix
          */
         TOperator& rRHSMatrix(const int idx);

         /**
          * @brief Set time dependent part of LHS matrix
          *
          * @param idx Index of the matrix
          */
         TOperator& rTMatrix(const int idx);

         /**
          * @brief Add RHS and solution data storage
          * 
          * @param rows Number of rows of matrix
          * @param cols Number of columns required
          */
         virtual void addStorage(const int rows, const int cols);
         
      protected:
         /**
          * @brief Complex RHS operators of the timestepped equations
          */
         std::vector<TOperator>   mRHSMatrix;

         /**
          * @brief Time dependent part of LHS matrix
          */
         std::vector<TOperator>   mTMatrix;

         /**
          * @brief Storage for old nonlinear RHS
          */
         std::vector<TData>  mRHSOld;

      private:
   };

   /// Typedef for a real data real matrix timestepper
   typedef SparseTimestepper<SparseMatrix,Matrix> SparseRTimestepper;
   /// Typedef for a shared real data real matrix timestepper
   typedef SharedPtrMacro<SparseRTimestepper> SharedSparseRTimestepper;
   /// Typedef for a complex data real matrix timestepper
   typedef SparseTimestepper<SparseMatrix,DecoupledZMatrix> SparseRZTimestepper;
   /// Typedef for a shared real data real matrix timestepper
   typedef SharedPtrMacro<SparseRZTimestepper> SharedSparseRZTimestepper;
   /// Typedef for a complex data complex matrix timestepper
   typedef SparseTimestepper<SparseMatrixZ,MatrixZ> SparseZTimestepper;
   /// Typedef for a shared real data real matrix timestepper
   typedef SharedPtrMacro<SparseZTimestepper> SharedSparseZTimestepper;

   template <typename TOperator,typename TData> SparseTimestepper<TOperator,TData>::SparseTimestepper(const int start)
      : Solver::SparseLinearSolver<TOperator,TData>(start)
   {
   }

   template <typename TOperator,typename TData> SparseTimestepper<TOperator,TData>::~SparseTimestepper()
   {
   }

   template <typename TOperator,typename TData> void SparseTimestepper<TOperator,TData>::computeRHS(const int step)
   {
      int start = step*this->nSystem();

      TData tmp;

      if(TimeSchemeType::rhsNN(step) == 0.0)
      {
         for(size_t i = this->mZeroIdx; i < this->mRHSData.size(); i++)
         {
            internal::computeRHSNoMemory<TOperator,TData>(step, this->mRHSData.at(i), this->mRHSMatrix.at(i+start), this->mSolution.at(i), this->mRHSOld.at(i));
         }
      } else
      {
         for(size_t i = this->mZeroIdx; i < this->mRHSData.size(); i++)
         {
            internal::computeRHSMemory<TOperator,TData>(step, this->mRHSData.at(i), this->mRHSMatrix.at(i+start), this->mSolution.at(i), this->mRHSOld.at(i), tmp);
         }
      }
   }

   template <typename TOperator,typename TData> void  SparseTimestepper<TOperator,TData>::updateTimeMatrix(const MHDFloat lhsCoeff, const MHDFloat rhsCoeff, const int step)
   {
      // Set start offset
      int start = step*this->nSystem();

      // Loop over matrices within same step
      for(int i = 0; i < this->nSystem(); ++i)
      {
         // Get the number of nonzero elements in time dependence
         size_t nnz = this->mTMatrix.at(i).nonZeros();

         // Update LHS and RHS matrices
         size_t lhsJ = 0;
         size_t rhsJ = 0;
         for (size_t k=0; k< static_cast<size_t>(this->mTMatrix.at(i).outerSize()); ++k)
         {
            typename TOperator::InnerIterator lhsIt(this->mLHSMatrix.at(start+i),lhsJ);
            typename TOperator::InnerIterator rhsIt(this->mRHSMatrix.at(start+i),rhsJ);
            for(typename TOperator::InnerIterator timeIt(this->mTMatrix.at(i),k); timeIt; ++timeIt)
            {
               // Only keep going if nonzero elements are left
               if(nnz > 0)
               {
                  // Update LHS matrix
                  if(timeIt.col() == lhsIt.col())
                  {
                     if(timeIt.row() == lhsIt.row())
                     {
                        // Update values
                        lhsIt.valueRef() += lhsCoeff*timeIt.value();

                        // Update LHS iterators and counters
                        ++lhsIt;
                        if(!lhsIt)
                        {
                           lhsJ++;
                        }
                     }
                  }

                  // Update LHS matrix
                  if(timeIt.col() == rhsIt.col())
                  {
                     if(timeIt.row() == rhsIt.row())
                     {
                        // Update values
                        rhsIt.valueRef() += rhsCoeff*timeIt.value();

                        // Update RHS iterators and counters
                        ++rhsIt;
                        if(!rhsIt)
                        {
                           rhsJ++;
                        }
                     }
                  }

                  // Update nonzero counter
                  nnz--;
               } else
               {
                  break;
               }
            }
         }

         // Safety assert to make sure all values have been updated
         assert(nnz == 0);
      }
   }

   template <typename TOperator,typename TData> void SparseTimestepper<TOperator,TData>::initMatrices(const int n)
   {
      // Initialise base matrices
      Solver::SparseLinearSolver<TOperator,TData>::initMatrices(n);

      // Do not reinitialise if work already done by other field
      if(this->mRHSMatrix.size() == 0)
      {
         // Reserve space for the RHS matrices
         this->mRHSMatrix.reserve(n);

         // Initialise storage
         for(int i = 0; i < n; ++i)
         {
            // Create storage for LHS matrices
            this->mRHSMatrix.push_back(TOperator());
         }
      }
   }

   template <typename TOperator,typename TData> TOperator& SparseTimestepper<TOperator,TData>::rRHSMatrix(const int idx)
   {
      return this->mRHSMatrix.at(idx);
   }

   template <typename TOperator,typename TData> TOperator& SparseTimestepper<TOperator,TData>::rTMatrix(const int idx)
   {
      return this->mTMatrix.at(idx);
   }

   template <typename TOperator,typename TData> void SparseTimestepper<TOperator,TData>::addStorage(const int rows, const int cols)
   {
      // Assert for non zero rows and columns
      assert(rows > 0);
      assert(cols > 0);

      Solver::SparseLinearSolver<TOperator,TData>::addStorage(rows,cols);

      // Add storage for old RHS
      this->mRHSOld.push_back(TData(rows,cols));
      this->mRHSOld.back().setZero();

      // Add storage for the time matrix
      this->mTMatrix.push_back(TOperator());
   }

   namespace internal
   {
      template <typename TOperator,typename TData> inline void computeRHSNoMemory(const int step, TData& rRHS, const TOperator& mat, const TData& sol, TData& rOld)
      {
         rOld = rRHS;
         rRHS = mat*sol + TimeSchemeType::rhsN(step)*rRHS;
      }

      template <> inline void computeRHSNoMemory<SparseMatrix,DecoupledZMatrix>(const int step, DecoupledZMatrix& rRHS, const SparseMatrix& mat, const DecoupledZMatrix& sol, DecoupledZMatrix& rOld)
      {
         rOld.real() = rRHS.real();
         rRHS.real() = mat*sol.real() + TimeSchemeType::rhsN(step)*rRHS.real();

         rOld.imag() = rRHS.imag();
         rRHS.imag() = mat*sol.imag() + TimeSchemeType::rhsN(step)*rRHS.imag();
      }

      template <typename TOperator,typename TData> inline void computeRHSMemory(const int step, TData& rRHS, const TOperator& mat, const TData& sol, TData& rOld, TData& rTmp)
      {
         rTmp = rRHS;
         rRHS = mat*sol + TimeSchemeType::rhsN(step)*rRHS + TimeSchemeType::rhsNN(step)*rOld;
         rOld = rTmp;
      }

      template <> inline void computeRHSMemory<SparseMatrix,DecoupledZMatrix>(const int step, DecoupledZMatrix& rRHS, const SparseMatrix& mat, const DecoupledZMatrix& sol, DecoupledZMatrix& rOld, DecoupledZMatrix& rTmp)
      {
         // rTmp is only a temporary variable, using only real reduces memory usage
         rTmp.real() = rRHS.real();
         rRHS.real() = mat*sol.real() + TimeSchemeType::rhsN(step)*rRHS.real() + TimeSchemeType::rhsNN(step)*rOld.real();
         rOld.real() = rTmp.real();

         rTmp.real() = rRHS.imag();
         rRHS.imag() = mat*sol.imag() + TimeSchemeType::rhsN(step)*rRHS.imag() + TimeSchemeType::rhsNN(step)*rOld.imag();
         rOld.imag() = rTmp.real();
      }
   }
}
}

#endif // SPARSETIMESTEPPER_HPP
