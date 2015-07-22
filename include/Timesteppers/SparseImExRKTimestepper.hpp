/** 
 * @file SparseImExRKTimestepper.hpp
 * @brief Implementation of a templated (coupled) equation timestepper for Implicit-Explicit Runge-Kutta schemes.
 *
 * The implementation is based on Cavaglieri & Bewley, "Low-storage implicit/explicit Runge-Kutta schemes for the simulation of stiff high-dimensional ODE systems", JCP, 2015
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef SPARSEIMEXRKTIMESTEPPER_HPP
#define SPARSEIMEXRKTIMESTEPPER_HPP

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
      /**
       * @brief Compute z = y
       */
      template <typename TOperator,typename TData> void computeSet(TData& z, const TData& y);

      /**
       * @brief Compute z = y
       */
      void computeSet(DecoupledZMatrix& z, const DecouplexZMatrix& y);

      /**
       * @brief Compute z = A*y
       */
      template <typename TOperator,typename TData> void computeMV(TData& z, const SparseMatrix& mat, const TData& y);

      /**
       * @brief Compute z = A*y
       */
      void computeMV(DecoupledZMatrix& z, const SparseMatrix& mat, const DecouplexZMatrix& y);

      /**
       * @brief Compute z = a*x + b*y + z
       */
      template <typename TOperator,typename TData> void computeAXPBYPZ(TData& z, const MHDFloat a, const TData& x, const MHDFloat b, const TData& y);

      /**
       * @brief Compute z = a*x + b*y + z
       */
      void computeAXPBYPZ(DecoupledZMatrix& z, const MHDFloat a, const DecoupledZMatrix& x, const MHDFloat b, const DecouplexZMatrix& y);

      /**
       * @brief Compute y = x + a*z
       */
      template <typename TOperator,typename TData> void computeXPAY(TData& y, const TData& x, const MHDFloat a);

      /**
       * @brief Compute y = x + a*y
       */
      void computeXPAY(DecoupledZMatrix& y, const DecoupledZMatrix& x, const MHDFloat a);

      /**
       * @brief Compute z = x + a*y + b*z
       */
      template <typename TOperator,typename TData> void computeXPAYBZ(TData& z, const TData& x, const MHDFloat a, const TData& y, const MHDFloat b);

      /**
       * @brief Compute z = x + a*y + b*z
       */
      void computeXPAYBZ(DecoupledZMatrix& z, const DecoupledZMatrix& x, const MHDFloat a, const DecoupledZMatrix& y, const MHDFloat b);
   }

   /**
    * @brief Implementation of a templated (coupled) equation timestepper for Implicit-Explicit Runge-Kutta schemes
    */
   template <typename TOperator,typename TData> class SparseImExRKTimestepper: public Solver::SparseLinearSolver<TOperator,TData>
   {
      public:
         /**
          * @brief Constructor
          *
          * @param start   Starting index (for example without m=0)
          * @param time    Solver timing with respect to timestepping
          */
         SparseImExRKTimestepper(const int start, const SolveTiming::Id time);

         /**
          * @brief Destructor
          */
         virtual ~SparseImExRKTimestepper();

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
          * @brief Prepare fields for implicit solve
          *
          * @param step    Current substep
          */
         bool preSolve(const int step);

         /**
          * @brief Work on fields after implicit solve
          *
          * @param step    Current substep
          */
         void postSolve(const int step);

         /**
          * @brief Update the LHS matrix with new timedependence
          *
          * @param step          Timestep scheme substep
          */
         void updateTimeMatrix(const int step, const MHDFloat dt);

         /**
          * @brief Set RHS matrix at t_n
          *
          * @param idx Index of the matrix
          */
         TOperator& rRHSMatrix(const int idx);

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
         void buildOperators(const int idx, const TOperator& opA, const TOperator& opB, const TOperator& opC, const MHDFloat dt, const int size);
         
      protected:
         /**
          * @brief Current timestep
          */
         MHDFloat mDt;

         /**
          * @brief RHS operator
          */
         std::vector<TOperator>   mRHSMatrix;

         /**
          * @brief Storage for implicit solution piece
          */
         std::vector<TData>  mImSolution;

         /**
          * @brief Storage for explicit solution piece
          */
         std::vector<TData>  mExSolution;

         /**
          * @brief Storage for intermediate solution
          */
         std::vector<TData>  mIntSolution;

         /**
          * @brief Storage for lower order embedded scheme solution
          */
         std::vector<TData>  mErrSolution;

      private:
   };

   template <typename TOperator,typename TData> SparseImExRKTimestepper<TOperator,TData>::SparseImExRKTimestepper(const int start, const SolveTiming::Id time)
      : Solver::SparseLinearSolver<TOperator,TData>(start, time), mDt(-1.0)
   {
   }

   template <typename TOperator,typename TData> SparseImExRKTimestepper<TOperator,TData>::~SparseImExRKTimestepper()
   {
   }

   template <typename TOperator,typename TData> bool SparseImExRKTimestepper<TOperator,TData>::preSolve(const int step)
   {
      if(step > 0)
      {
         // Update intermediate solution
         MHDFloat bIm = TimeSchemeSelector::bIm(step)*this->mDt;
         MHDFloat bEx = TimeSchemeSelector::bEx(step)*this->mDt;
         for(size_t i = this->mZeroIdx; i < this->mRHSData.size(); i++)
         {
            internal::computeAXPBYPZ(this->mIntSolution.at(i), bIm, this->mImSolution.at(i), bEx, this->mExSolution.at(i));
         }

         // Embedded lower order scheme solution
         bIm = TimeSchemeSelector::bImErr(step)*this->mDt;
         bEx = TimeSchemeSelector::bExErr(step)*this->mDt;
         for(size_t i = this->mZeroIdx; i < this->mRHSData.size(); i++)
         {
            internal::computeAXPBYPZ(this->mErrSolution.at(i), bIm, this->mImSolution.at(i), bEx, this->mExSolution.at(i));
         }
      }

      // Update explicit term
      for(size_t i = this->mZeroIdx; i < this->mRHSData.size(); i++)
      {
         internal::computeSet(this->mExSolution.at(i), this->mRHSData.at(i));
      }

      // First step
      if(step == 0)
      {
         // Build RHS for implicit term
         for(size_t i = this->mZeroIdx; i < this->mRHSData.size(); i++)
         {
            internal::computeMV(this->mRHSData.at(i), this->mRHSMatrix.at(i), this->mIntSolution.at(i));
         }

         return true;
      
      // Last step has no implicit solve
      } else if(step == TimeSchemeSelector::STEPS - 1)
      {
         for(size_t i = this->mZeroIdx; i < this->mRHSData.size(); i++)
         {
            internal::computeSet(this->mSolution.at(i), this->mIntSolution.at(i));
         }

         return false;

      } else
      {
         // Build RHS for implicit term
         MHDFloat aIm = (TimeSchemeSelector::aIm(step, step-1) - TimeSchemeSelector::bIm(step-1))*this->mDt;
         MHDFloat aEx = (TimeSchemeSelector::aEx(step, step-1) - TimeSchemeSelector::bEx(step-1))*this->mDt;
         for(size_t i = this->mZeroIdx; i < this->mRHSData.size(); i++)
         {
            internal::computeXPAYPBZ(this->mExSolution.at(i), this->mIntSolution.at(i), aIm, this->mImSolution.at(i), aEx);
            internal::computeMV(this->mRHSData.at(i), this->mRHSMatrix.at(i), this->mExSolution.at(i));
         }

         return true
      }
   }

   template <typename TOperator,typename TData> void SparseImExRKTimestepper<TOperator,TData>::postSolve(const int step)
   {
      // Update implicit term
      for(size_t i = this->mZeroIdx; i < this->mRHSData.size(); i++)
      {
         this->mImSolution.at(i) = this->mSolution.at(i);
      }

      if(step == 0)
      {
         // Prepare solution for new nonlinear term
         MHDFloat aIm = TimeSchemeSelector::aIm(step, step)*this->mDt;
         for(size_t i = this->mZeroIdx; i < this->mRHSData.size(); i++)
         {
            internal::computeXPAY(this->mSolution.at(i), this->mIntSolution.at(i), aIm);
         }
      } else
      {
         // Prepare solution for new nonlinear term
         MHDFloat aIm = TimeSchemeSelector::aIm(step, step)*this->mDt;
         for(size_t i = this->mZeroIdx; i < this->mRHSData.size(); i++)
         {
            internal::computeXPAY(this->mSolution.at(i), this->mExSolution.at(i), aIm);
         }
      }
   }

   template <typename TOperator,typename TData> void  SparseImExRKTimestepper<TOperator,TData>::updateTimeMatrix(const int step, const MHDFloat dt)
   {
      // Update stored timestep
      this->mDt = dt;

      // Set start offset
      int start = step*this->nSystem();

      // Loop over matrices within same step
      for(int i = 0; i < this->nSystem(); ++i)
      {
         // Get the number of nonzero elements in time dependence
         size_t nnz = this->mRHSMatrix.at(i).nonZeros();

         // Update LHS and RHS matrices
         size_t lhsJ = 0;
         size_t rhsJ = 0;
         for (size_t k=0; k< static_cast<size_t>(this->mRHSMatrix.at(i).outerSize()); ++k)
         {
            typename TOperator::InnerIterator lhsIt(this->mLHSMatrix.at(start+i),lhsJ);
            typename TOperator::InnerIterator rhsIt(this->mRHSMatrix.at(start+i),rhsJ);
            for(typename TOperator::InnerIterator timeIt(this->mRHSMatrix.at(i),k); timeIt; ++timeIt)
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

                  // Update RHS matrix at t_n
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

   template <typename TOperator,typename TData> void SparseImExRKTimestepper<TOperator,TData>::initMatrices(const int n)
   {
      // Initialise base matrices
      Solver::SparseLinearSolver<TOperator,TData>::initMatrices(TimeSchemeSelector::STEPS*n);

      // Do not reinitialise if work already done by other field
      if(this->mRHSMatrix.size() == 0)
      {
         // Reserve space for the RHS matrices
         this->mRHSMatrix.reserve(n);

         // Initialise storage for RHS matrices
         for(int i = 0; i < n; ++i)
         {
            // Create storage for LHS matrices
            this->mRHSMatrix.push_back(TOperator());
         }
      }
   }

   template <typename TOperator,typename TData> void SparseImExRKTimestepper<TOperator,TData>::initSolutions()
   {
      for(size_t i = this->mZeroIdx; i < this->mRHSData.size(); i++)
      {
         internal::computeSet(this->mIntSolution.at(i), this->mSolution.at(i));
      }
   }

   template <typename TOperator,typename TData> TOperator& SparseImExRKTimestepper<TOperator,TData>::rRHSMatrix(const int idx)
   {
      return this->mRHSMatrix.at(idx);
   }

   template <typename TOperator,typename TData> void SparseImExRKTimestepper<TOperator,TData>::addStorage(const int rows, const int cols)
   {
      // Assert for non zero rows and columns
      assert(rows > 0);
      assert(cols > 0);

      Solver::SparseLinearSolver<TOperator,TData>::addStorage(rows,cols);

      // Add RK storage
      this->mImSolution.push_back(TData(rows,cols));
      this->mExSolution.push_back(TData(rows,cols));
      this->mIntSolution.push_back(TData(rows,cols));
      this->mErrSolution.push_back(TData(rows,cols));
   }

   template <typename TOperator,typename TData> void SparseImExRKTimestepper<TOperator,TData>::buildOperators(const int idx, const TOperator& opA, const TOperator& opB, const TOperator& opC, const MHDFloat dt, const int size)
   {
      // Update timestep
      this->mDt = dt;

      // Resize matrices if necessary
      this->rLHSMatrix(idx).resize(size, size);
      this->rRHSMatrix(idx).resize(size, size);

      this->rRHSMatrix(idx) = opA;
      for(int i = 0; i < TimeSchemeSelector::STEPS; ++i)
      {
         this->rLHSMatrix(idx + i*this->nSystem()) = opB - TimeSchemeSelector::aIm(i,i)*this->mDt*opA + opC;
      }
   }

   namespace internal
   {
      template <typename TOperator,typename TData> inline void computeAXPBYPZ(TData& y, const TData& x)
      {
         z = y;
      }

      inline void computeAXPBYPZ(DecoupledZMatrix& z, const DecoupledZMatrix& x)
      {
         y.real() = x.real();

         y.imag() = x.imag();
      }

      template <typename TOperator,typename TData> inline void computeXPAY(TData& y, const TData& x, const MHDFloat a)
      {
         y = x + a*y;
      }

      inline void computeXPAY(DecoupledZMatrix& y, const DecoupledZMatrix& x, const MHDFloat a)
      {
         y.real() = x.real() + a*y.real();

         y.imag() = x.imag() + a*y.imag();
      }

      template <typename TOperator,typename TData> inline void computeXPAYPBZ(TData& z, const TData& x, const MHDFloat a, const TData& y, const MHDFloat b)
      {
        z = x + a*y + b*z;
      }

      inline void computeXPAYPBZ(DecoupledZMatrix& z, const DecoupledZMatrix& x, const MHDFloat a, const DecoupledZMatrix& y, const MHDFloat b)
      {
         z.real() = x.real() + a*y.real() + b*z.real();

         z.imag() = x.imag() + a*y.imag() + b*z.imag();
      }

      template <typename TOperator,typename TData> inline void computeAXPBYPZ(TData& z, const MHDFloat a, const TData& x, const MHDFloat b, const TData& y)
      {
         z += a*x + b*y;
      }

      inline void computeAXPBYPZ(DecoupledZMatrix& z, const MHDFloat a, const DecoupledZMatrix& x, const MHDFloat b, const DecoupledZMatrix& y)
      {
         z.real() += a*x.real() + b*y.real();

         z.imag() += a*x.imag() + b*y.imag();
      }

      template <typename TOperator,typename TData> inline void computeMV(TData& y, const SparseMatrix& A, const TData& x)
      {
         y = A*x;
      }

      inline void computeMV(DecoupledZMatrix& y, const SparseMatrix& A, const DecoupledZMatrix& x)
      {
         y.real() = A*x.real();

         y.imag() = A*x.imag();
      }
   }
}
}

#endif // SPARSEIMEXRKTIMESTEPPER_HPP
