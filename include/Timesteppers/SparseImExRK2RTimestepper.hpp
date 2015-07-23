/** 
 * @file SparseImExRK2RTimestepper.hpp
 * @brief Implementation of a templated (coupled) equation timestepper for Implicit-Explicit Runge-Kutta (2R) schemes.
 *
 * The implementation is based on Cavaglieri & Bewley, "Low-storage implicit/explicit Runge-Kutta schemes for the simulation of stiff high-dimensional ODE systems", JCP, 2015
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef SPARSEIMEXRK2RTIMESTEPPER_HPP
#define SPARSEIMEXRK2RTIMESTEPPER_HPP

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
      template <typename TData> void computeSet(TData& z, const TData& y);

      /**
       * @brief Compute z = y
       */
      void computeSet(DecoupledZMatrix& z, const DecoupledZMatrix& y);

      /**
       * @brief Compute z = A*y
       */
      template <typename TOperator,typename TData> void computeMV(TData& z, const TOperator& mat, const TData& y);

      /**
       * @brief Compute z = A*y
       */
      void computeMV(DecoupledZMatrix& z, const SparseMatrix& mat, const DecoupledZMatrix& y);

      /**
       * @brief Compute z = a*x + b*y + z
       */
      template <typename TData> void computeAXPBYPZ(TData& z, const MHDFloat a, const TData& x, const MHDFloat b, const TData& y);

      /**
       * @brief Compute z = a*x + b*y + z
       */
      void computeAXPBYPZ(DecoupledZMatrix& z, const MHDFloat a, const DecoupledZMatrix& x, const MHDFloat b, const DecoupledZMatrix& y);

      /**
       * @brief Compute y = x + a*z
       */
      template <typename TData> void computeXPAY(TData& y, const TData& x, const MHDFloat a);

      /**
       * @brief Compute y = x + a*y
       */
      void computeXPAY(DecoupledZMatrix& y, const DecoupledZMatrix& x, const MHDFloat a);

      /**
       * @brief Compute z = x + a*y + b*z
       */
      template <typename TData> void computeXPAYPBZ(TData& z, const TData& x, const MHDFloat a, const TData& y, const MHDFloat b);

      /**
       * @brief Compute z = x + a*y + b*z
       */
      void computeXPAYPBZ(DecoupledZMatrix& z, const DecoupledZMatrix& x, const MHDFloat a, const DecoupledZMatrix& y, const MHDFloat b);
   }

   /**
    * @brief Implementation of a templated (coupled) equation timestepper for Implicit-Explicit Runge-Kutta (2R) schemes
    */
   template <typename TOperator,typename TData> class SparseImExRK2RTimestepper: public Solver::SparseLinearSolver<TOperator,TData>
   {
      public:
         /**
          * @brief Constructor
          *
          * @param start   Starting index (for example without m=0)
          * @param time    Solver timing with respect to timestepping
          */
         SparseImExRK2RTimestepper(const int start, const SolveTiming::Id time);

         /**
          * @brief Destructor
          */
         virtual ~SparseImExRK2RTimestepper();

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
          */
         bool preSolve();

         /**
          * @brief Work on fields after implicit solve
          *
          * @param step    Current substep
          */
         void postSolve();

         /**
          * @brief Update the LHS matrix with new timedependence
          */
         void updateTimeMatrix(const MHDFloat dt);

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
         void buildOperators(const int idx, const DecoupledZSparse& opA, const DecoupledZSparse& opB, const DecoupledZSparse& opC, const MHDFloat dt, const int size);

         /**
          * @brief Finished timestep?
          */
         bool finished();
         
      protected:
         /**
          * @brief Current substep
          */
         int mStep;

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

   template <typename TOperator,typename TData> SparseImExRK2RTimestepper<TOperator,TData>::SparseImExRK2RTimestepper(const int start, const SolveTiming::Id time)
      : Solver::SparseLinearSolver<TOperator,TData>(start, time), mStep(0), mDt(-1.0)
   {
   }

   template <typename TOperator,typename TData> SparseImExRK2RTimestepper<TOperator,TData>::~SparseImExRK2RTimestepper()
   {
   }

   template <typename TOperator,typename TData> bool SparseImExRK2RTimestepper<TOperator,TData>::finished()
   {
      this->mStep = (this->mStep + 1) % (TimeSchemeSelector::STEPS + 1);

      return (this->mStep == 0);
   }

   template <typename TOperator,typename TData> bool SparseImExRK2RTimestepper<TOperator,TData>::preSolve()
   {
      this->mId = TimeSchemeSelector::aIm(this->mStep, this->mStep);

      if(this->mStep > 0)
      {
         // Update intermediate solution
         MHDFloat bIm = TimeSchemeSelector::bIm(this->mStep-1)*this->mDt;
         MHDFloat bEx = TimeSchemeSelector::bEx(this->mStep-1)*this->mDt;
         for(size_t i = this->mZeroIdx; i < this->mRHSData.size(); i++)
         {
            internal::computeAXPBYPZ(this->mIntSolution.at(i), bIm, this->mImSolution.at(i), bEx, this->mExSolution.at(i));
         }

         // Embedded lower order scheme solution
         if(TimeSchemeSelector::HAS_EMBEDDED)
         {
            bIm = TimeSchemeSelector::bImErr(this->mStep-1)*this->mDt;
            bEx = TimeSchemeSelector::bExErr(this->mStep-1)*this->mDt;
            for(size_t i = this->mZeroIdx; i < this->mRHSData.size(); i++)
            {
               internal::computeAXPBYPZ(this->mErrSolution.at(i), bIm, this->mImSolution.at(i), bEx, this->mExSolution.at(i));
            }
         }
      }

      // Update explicit term
      for(size_t i = this->mZeroIdx; i < this->mRHSData.size(); i++)
      {
         internal::computeSet(this->mExSolution.at(i), this->mRHSData.at(i));
      }

      // First step
      if(this->mStep == 0)
      {
         // Build RHS for implicit term
         for(size_t i = this->mZeroIdx; i < this->mRHSData.size(); i++)
         {
            internal::computeMV(this->mRHSData.at(i), this->mRHSMatrix.at(i), this->mIntSolution.at(i));
         }

         return true;
      
      // Last step has no implicit solve
      } else if(this->mStep == TimeSchemeSelector::STEPS)
      {
         for(size_t i = this->mZeroIdx; i < this->mRHSData.size(); i++)
         {
            internal::computeSet(this->mSolution.at(i), this->mIntSolution.at(i));
         }

         return false;

      } else
      {
         // Build RHS for implicit term
         MHDFloat aIm = (TimeSchemeSelector::aIm(this->mStep, this->mStep-1) - TimeSchemeSelector::bIm(this->mStep-1))*this->mDt;
         MHDFloat aEx = (TimeSchemeSelector::aEx(this->mStep, this->mStep-1) - TimeSchemeSelector::bEx(this->mStep-1))*this->mDt;
         for(size_t i = this->mZeroIdx; i < this->mRHSData.size(); i++)
         {
            internal::computeXPAYPBZ(this->mExSolution.at(i), this->mIntSolution.at(i), aIm, this->mImSolution.at(i), aEx);
            internal::computeMV(this->mRHSData.at(i), this->mRHSMatrix.at(i), this->mExSolution.at(i));
         }

         return true;
      }
   }

   template <typename TOperator,typename TData> void SparseImExRK2RTimestepper<TOperator,TData>::postSolve()
   {
      // Update implicit term
      for(size_t i = this->mZeroIdx; i < this->mRHSData.size(); i++)
      {
         this->mImSolution.at(i) = this->mSolution.at(i);
      }

      if(this->mStep == 0)
      {
         // Prepare solution for new nonlinear term
         MHDFloat aIm = TimeSchemeSelector::aIm(this->mStep, this->mStep)*this->mDt;
         for(size_t i = this->mZeroIdx; i < this->mRHSData.size(); i++)
         {
            internal::computeXPAY(this->mSolution.at(i), this->mIntSolution.at(i), aIm);
         }
      } else
      {
         // Prepare solution for new nonlinear term
         MHDFloat aIm = TimeSchemeSelector::aIm(this->mStep, this->mStep)*this->mDt;
         for(size_t i = this->mZeroIdx; i < this->mRHSData.size(); i++)
         {
            internal::computeXPAY(this->mSolution.at(i), this->mExSolution.at(i), aIm);
         }
      }
   }

   template <typename TOperator,typename TData> void  SparseImExRK2RTimestepper<TOperator,TData>::updateTimeMatrix(const MHDFloat dt)
   {
      // Update stored timestep
      MHDFloat oldDt = this->mDt;
      this->mDt = dt;

      for(int step = 0; step < TimeSchemeSelector::STEPS; ++step)
      {
         // Loop over matrices within same step
         for(int i = 0; i < this->nSystem(); ++i)
         {
            // Get the number of nonzero elements in time dependence
            size_t nnz = this->mRHSMatrix.at(i).nonZeros();

            // Update LHS and RHS matrices
            size_t lhsJ = 0;
            for (size_t k=0; k< static_cast<size_t>(this->mRHSMatrix.at(i).outerSize()); ++k)
            {
               typename TOperator::InnerIterator lhsIt(this->rLHSMatrix(TimeSchemeSelector::aIm(step,step), i),lhsJ);
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
                           lhsIt.valueRef() += TimeSchemeSelector::aIm(step,step)*(oldDt - this->mDt)*timeIt.value();

                           // Update LHS iterators and counters
                           ++lhsIt;
                           if(!lhsIt)
                           {
                              lhsJ++;
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
   }

   template <typename TOperator,typename TData> void SparseImExRK2RTimestepper<TOperator,TData>::initMatrices(const int n)
   {
      // Initialise base matrices
      for(int i = 0; i < TimeSchemeSelector::STEPS; i++)
      {
         Solver::SparseLinearSolver<TOperator,TData>::initMatrices(TimeSchemeSelector::aIm(i,i), n);
      }

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

   template <typename TOperator,typename TData> void SparseImExRK2RTimestepper<TOperator,TData>::initSolutions()
   {
      for(size_t i = this->mZeroIdx; i < this->mRHSData.size(); i++)
      {
         internal::computeSet(this->mIntSolution.at(i), this->mSolution.at(i));
      }
   }

   template <typename TOperator,typename TData> TOperator& SparseImExRK2RTimestepper<TOperator,TData>::rRHSMatrix(const int idx)
   {
      return this->mRHSMatrix.at(idx);
   }

   template <typename TOperator,typename TData> void SparseImExRK2RTimestepper<TOperator,TData>::addStorage(const int rows, const int cols)
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

   template <typename TOperator,typename TData> void SparseImExRK2RTimestepper<TOperator,TData>::buildOperators(const int idx, const DecoupledZSparse& opA, const DecoupledZSparse& opB, const DecoupledZSparse& opC, const MHDFloat dt, const int size)
   {
      // Update timestep
      this->mDt = dt;

      // Set explicit matrix
      this->rRHSMatrix(idx).resize(size, size);
      Solver::internal::addOperators(this->rRHSMatrix(idx), 1.0, opA);

      // Set implicit matrix
      for(int i = 0; i < TimeSchemeSelector::STEPS; ++i)
      {
         this->rLHSMatrix(TimeSchemeSelector::aIm(i,i), idx).resize(size, size);
         Solver::internal::addOperators(this->rLHSMatrix(TimeSchemeSelector::aIm(i,i), idx), 1.0, opB);
         Solver::internal::addOperators(this->rLHSMatrix(TimeSchemeSelector::aIm(i,i), idx), -TimeSchemeSelector::aIm(i,i)*this->mDt, opA);
         Solver::internal::addOperators(this->rLHSMatrix(TimeSchemeSelector::aIm(i,i), idx), 1.0, opC);
      }
   }

   namespace internal
   {
      template <typename TData> inline void computeSet(TData& y, const TData& x)
      {
         y = x;
      }

      inline void computeSet(DecoupledZMatrix& y, const DecoupledZMatrix& x)
      {
         y.real() = x.real();

         y.imag() = x.imag();
      }

      template <typename TData> inline void computeXPAY(TData& y, const TData& x, const MHDFloat a)
      {
         y = x + a*y;
      }

      inline void computeXPAY(DecoupledZMatrix& y, const DecoupledZMatrix& x, const MHDFloat a)
      {
         y.real() = x.real() + a*y.real();

         y.imag() = x.imag() + a*y.imag();
      }

      template <typename TData> inline void computeXPAYPBZ(TData& z, const TData& x, const MHDFloat a, const TData& y, const MHDFloat b)
      {
        z = x + a*y + b*z;
      }

      inline void computeXPAYPBZ(DecoupledZMatrix& z, const DecoupledZMatrix& x, const MHDFloat a, const DecoupledZMatrix& y, const MHDFloat b)
      {
         z.real() = x.real() + a*y.real() + b*z.real();

         z.imag() = x.imag() + a*y.imag() + b*z.imag();
      }

      template <typename TData> inline void computeAXPBYPZ(TData& z, const MHDFloat a, const TData& x, const MHDFloat b, const TData& y)
      {
         z += a*x + b*y;
      }

      inline void computeAXPBYPZ(DecoupledZMatrix& z, const MHDFloat a, const DecoupledZMatrix& x, const MHDFloat b, const DecoupledZMatrix& y)
      {
         z.real() += a*x.real() + b*y.real();

         z.imag() += a*x.imag() + b*y.imag();
      }

      template <typename TOperator,typename TData> inline void computeMV(TData& y, const TOperator& A, const TData& x)
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

#endif // SPARSEIMEXRK2RTIMESTEPPER_HPP
