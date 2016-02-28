/** 
 * @file SparseImExRK3RTimestepper.hpp
 * @brief Implementation of a templated (coupled) equation timestepper for Implicit-Explicit Runge-Kutta (3R) schemes.
 *
 * The implementation is based on Cavaglieri & Bewley, "Low-storage implicit/explicit Runge-Kutta schemes for the simulation of stiff high-dimensional ODE systems", JCP, 2015
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef SPARSEIMEXRK3RTIMESTEPPER_HPP
#define SPARSEIMEXRK3RTIMESTEPPER_HPP

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
       * @brief Compute z = a*y
       */
      template <typename TData> void computeSet(TData& z, const MHDFloat a, const TData& y);

      /**
       * @brief Compute z = a*y
       */
      void computeSet(DecoupledZMatrix& z, const MHDFloat a, const DecoupledZMatrix& y);

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
       * @brief Compute y = a*M*x + y
       */
      template <typename TData> void computeAMXPY(TData& y, const SparseMatrix& mat, const MHDFloat a, const TData& x);

      /**
       * @brief Compute y = a*M*x + y
       */
      void computeAMXPY(DecoupledZMatrix& y, const SparseMatrix& mat, const MHDFloat a, const DecoupledZMatrix& x);

      /**
       * @brief Compute z = a*M*x + y + b*z
       */
      template <typename TData> void computeAMXPYPBZ(TData& z, const SparseMatrix& mat, const MHDFloat a, const TData& x, const TData& y, const MHDFloat b);

      /**
       * @brief Compute z = a*M*x + y + b*z
       */
      void computeAMXPYPBZ(DecoupledZMatrix& z, const SparseMatrix& mat, const MHDFloat a, const DecoupledZMatrix& x, const DecoupledZMatrix& y, const MHDFloat b);

      /**
       * @brief Compute z = a*M*x + b*y + z
       */
      template <typename TData> void computeAMXPBYPZ(TData& z, const SparseMatrix& mat, const MHDFloat a, const TData& x, const MHDFloat b, const TData& y);

      /**
       * @brief Compute z = a*M*x + b*y + z
       */
      void computeAMXPBYPZ(DecoupledZMatrix& z, const SparseMatrix& mat, const MHDFloat a, const DecoupledZMatrix& x, const MHDFloat b, const DecoupledZMatrix& y);

      /**
       * @brief Compute z = a*M*x + b*y + M*z
       */
      template <typename TData> void computeAMXPBYPMZ(TData& z, const SparseMatrix& mat, const MHDFloat a, const TData& x, const MHDFloat b, const TData& y);

      /**
       * @brief Compute z = a*M*x + b*y + M*z
       */
      void computeAMXPBYPMZ(DecoupledZMatrix& z, const SparseMatrix& mat, const MHDFloat a, const DecoupledZMatrix& x, const MHDFloat b, const DecoupledZMatrix& y);

      /**
       * @brief Compute y = a*x + y
       */
      template <typename TData> void computeAXPY(TData& y, const MHDFloat a, const TData& x);

      /**
       * @brief Compute y = a*x + y
       */
      void computeAXPY(DecoupledZMatrix& y, const MHDFloat a, const DecoupledZMatrix& x);

      /**
       * @brief Compute y = x + a*y
       */
      template <typename TData> void computeXPAY(TData& y, const TData& x, const MHDFloat a);

      /**
       * @brief Compute y = x + a*y
       */
      void computeXPAY(DecoupledZMatrix& y, const DecoupledZMatrix& x, const MHDFloat a);

      /**
       * @brief Compute z = M*y
       */
      template <typename TOperator,typename TData> void computeMY(TData& z, const TOperator& mat, const TData& y);

      /**
       * @brief Compute z = M*y
       */
      void computeMY(DecoupledZMatrix& z, const SparseMatrix& mat, const DecoupledZMatrix& y);

      /**
       * @brief Compute error
       */
      template <typename TData> void computeError(MHDFloat& err, const TData& x, const TData& y);

      /**
       * @brief Compute error
       */
      void computeError(MHDFloat& err, const DecoupledZMatrix& x, const DecoupledZMatrix& y);

   }

   /**
    * @brief Implementation of a templated (coupled) equation timestepper for Implicit-Explicit Runge-Kutta (3R) schemes
    */
   template <typename TOperator,typename TData> class SparseImExRK3RTimestepper: public Solver::SparseLinearSolver<TOperator,TData>
   {
      public:
         /**
          * @brief Enum for registers
          */
         enum RegisterId {
            // 
            IMPLICIT_REGISTER = 0,
            // 
            EXPLICIT_REGISTER,
            // 
            TEMPORARY_REGISTER,
            // 
            SOLUTION_REGISTER,
            // 
            ERROR_REGISTER,
         };

         /**
          * @brief Constructor
          *
          * @param start   Starting index (for example without m=0)
          * @param time    Solver timing with respect to timestepping
          */
         SparseImExRK3RTimestepper(const int start, const SolveTiming::Id time);

         /**
          * @brief Destructor
          */
         virtual ~SparseImExRK3RTimestepper();

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
         virtual bool preSolve();

         /**
          * @brief Work on fields after implicit solve
          *
          * @param step    Current substep
          */
         virtual bool postSolve();

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
         MHDFloat error() const;

         /**
          * @brief Finished timestep?
          */
         bool finished();

         /**
          * @brief Get current timestep fraction
          */
         MHDFloat stepFraction() const;

      protected:
         /**
          * @brief Explicit calculation took place?
          */
         bool mHasExplicit;

         /**
          * @brief Current substep
          */
         int mStep;

         /**
          * @brief Current timestep
          */
         MHDFloat mDt;

         /**
          * @brief Timestep error
          */
         MHDFloat mError;

         /**
          * @brief ID of the register to use
          */
         RegisterId  mRegisterId;

         /**
          * @brief RHS operator
          */
         std::vector<TOperator>   mRHSMatrix;

         /**
          * @brief Mass matrix operator
          */
         std::vector<SparseMatrix>   mMassMatrix;

         /**
          * @brief Storage for implicit solution piece
          */
         std::vector<TData>  mImSolution;

         /**
          * @brief Storage for explicit solution piece
          */
         std::vector<TData>  mExSolution;

         /**
          * @brief Storage for temporary solution
          */
         std::vector<TData>  mTmpSolution;

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

   template <typename TOperator,typename TData> SparseImExRK3RTimestepper<TOperator,TData>::SparseImExRK3RTimestepper(const int start, const SolveTiming::Id time)
      : Solver::SparseLinearSolver<TOperator,TData>(start, time), mHasExplicit(true), mStep(0), mDt(-1.0), mError(-1.0), mRegisterId(IMPLICIT_REGISTER)
   {
   }

   template <typename TOperator,typename TData> SparseImExRK3RTimestepper<TOperator,TData>::~SparseImExRK3RTimestepper()
   {
   }

   template <typename TOperator,typename TData> MHDFloat SparseImExRK3RTimestepper<TOperator,TData>::error() const
   {
      return this->mError;
   }

   template <typename TOperator,typename TData> bool SparseImExRK3RTimestepper<TOperator,TData>::finished()
   {
      return (this->mStep == 0);
   }

   template <typename TOperator,typename TData> MHDFloat SparseImExRK3RTimestepper<TOperator,TData>::stepFraction() const
   {
      return TimeSchemeSelector::cEx(this->mStep);
   }

   template <typename TOperator,typename TData> bool SparseImExRK3RTimestepper<TOperator,TData>::preSolve()
   {
      if(this->mHasExplicit)
      {
         // Update explicit term with explicit (nonlinear) values
         for(size_t i = this->mZeroIdx; i < this->mRHSData.size(); i++)
         {
            internal::computeSet(this->mExSolution.at(i), this->mRHSData.at(i));
         }
      }

      if(this->mHasExplicit && this->mStep > 0)
      {
         // Update intermediate solution
         MHDFloat bIm = TimeSchemeSelector::bIm(this->mStep)*this->mDt;
         MHDFloat bEx = TimeSchemeSelector::bEx(this->mStep)*this->mDt;
         for(size_t i = this->mZeroIdx; i < this->mRHSData.size(); i++)
         {
            internal::computeAMXPBYPZ(this->mIntSolution.at(i), this->mMassMatrix.at(i), bIm, this->mImSolution.at(i), bEx, this->mExSolution.at(i));
         }

         // Embedded lower order scheme solution
         if(TimeSchemeSelector::HAS_EMBEDDED)
         {
            bIm = TimeSchemeSelector::bImErr(this->mStep)*this->mDt;
            bEx = TimeSchemeSelector::bExErr(this->mStep)*this->mDt;
            for(size_t i = this->mZeroIdx; i < this->mRHSData.size(); i++)
            {
               internal::computeAMXPBYPZ(this->mErrSolution.at(i), this->mMassMatrix.at(i), bIm, this->mImSolution.at(i), bEx, this->mExSolution.at(i));
            }
         }

         this->mStep += 1;
      }

      // First step
      if(this->mStep == 0)
      {
         // Reset error
         if(TimeSchemeSelector::HAS_EMBEDDED)
         {
            this->mError = 0.0;
         }

         // Build RHS for implicit term
         for(size_t i = this->mZeroIdx; i < this->mRHSData.size(); i++)
         {
            internal::computeMV(this->mRHSData.at(i), this->mRHSMatrix.at(i), this->mIntSolution.at(i));
         }

         // Initialise temporary storage
         for(size_t i = this->mZeroIdx; i < this->mRHSData.size(); i++)
         {
            internal::computeMV(this->mTmpSolution.at(i), this->mMassMatrix.at(i), this->mIntSolution.at(i));
         }

         // Set ID for solver
         this->mId = TimeSchemeSelector::aIm(this->mStep, this->mStep);

         this->mRegisterId = IMPLICIT_REGISTER;

         return true;
      
      // Last step has no implicit solve
      } else if(this->mStep == TimeSchemeSelector::STEPS)
      {
         if(this->mRegisterId == SOLUTION_REGISTER)
         {
            // Compute error estimate using embedded scheme
            if(TimeSchemeSelector::HAS_EMBEDDED)
            {
               for(size_t i = this->mZeroIdx; i < this->mRHSData.size(); i++)
               {
                  internal::computeSet(this->mRHSData.at(i), this->mErrSolution.at(i));
               }

               // Set mass matrix ID for solver
               this->mId = 0.0;

               // Set explicit store register for solution
               this->mRegisterId = ERROR_REGISTER; 
               
               return true;
            } else
            {
               for(size_t i = this->mZeroIdx; i < this->mRHSData.size(); i++)
               {
                  internal::computeSet(this->mSolution.at(i), this->mIntSolution.at(i));
               }

               // Explicit nonlinear term at next step
               this->mHasExplicit = true;

               // Reset step to 0
               this->mStep = 0;

               // Reset register ID
               this->mRegisterId = IMPLICIT_REGISTER; 

               return false;
            }
         } else if(this->mRegisterId == ERROR_REGISTER)
         {
            for(size_t i = this->mZeroIdx; i < this->mRHSData.size(); i++)
            {
               internal::computeSet(this->mSolution.at(i), this->mIntSolution.at(i));
            }

            for(size_t i = this->mZeroIdx; i < this->mRHSData.size(); i++)
            {
               internal::computeError(this->mError, this->mIntSolution.at(i), this->mErrSolution.at(i));

               internal::computeSet(this->mErrSolution.at(i), this->mIntSolution.at(i));
            }

            // Explicit nonlinear term at next step
            this->mHasExplicit = true;

            // Reset step to 0
            this->mStep = 0;

            // Reset register ID
            this->mRegisterId = IMPLICIT_REGISTER; 

            return false;
         } else
         {
            for(size_t i = this->mZeroIdx; i < this->mRHSData.size(); i++)
            {
               internal::computeSet(this->mRHSData.at(i), this->mIntSolution.at(i));
            }

            // Set mass matrix ID for solver
            this->mId = 0.0;

            // Set explicit store register for solution
            this->mRegisterId = SOLUTION_REGISTER;

            return true;
         }
      } else
      {
         if(this->mRegisterId == IMPLICIT_REGISTER)
         {
            MHDFloat aIm = 0.0;
            MHDFloat aEx = 0.0;

            // Build explicit term
            aEx = TimeSchemeSelector::aEx(this->mStep, this->mStep-1)*this->mDt;
            for(size_t i = this->mZeroIdx; i < this->mRHSData.size(); i++)
            {
               internal::computeXPAY(this->mExSolution.at(i), this->mTmpSolution.at(i), aEx);
            }

            if(this->mStep < TimeSchemeSelector::STEPS-1)
            {
               // Build RHS for implicit term
               aIm = (TimeSchemeSelector::aIm(this->mStep+1, this->mStep-1) - TimeSchemeSelector::bIm(this->mStep-1))*this->mDt;
               aEx = (TimeSchemeSelector::aEx(this->mStep+1, this->mStep-1) - TimeSchemeSelector::bEx(this->mStep-1))/TimeSchemeSelector::aEx(this->mStep, this->mStep-1);
               for(size_t i = this->mZeroIdx; i < this->mRHSData.size(); i++)
               {
                  internal::computeXPAY(this->mTmpSolution.at(i), this->mExSolution.at(i), -1.0);
                  internal::computeAMXPYPBZ(this->mTmpSolution.at(i), this->mMassMatrix.at(i), aIm, this->mImSolution.at(i), this->mIntSolution.at(i), aEx);
               }
            }

            // Build explicit term
            aIm = TimeSchemeSelector::aIm(this->mStep, this->mStep-1)*this->mDt;
            for(size_t i = this->mZeroIdx; i < this->mRHSData.size(); i++)
            {
               internal::computeAMXPY(this->mExSolution.at(i), this->mMassMatrix.at(i), aIm, this->mImSolution.at(i));
               internal::computeSet(this->mRHSData.at(i), this->mExSolution.at(i));
            }

            // Set mass matrix ID for solver
            this->mId = 0.0;

            // Set explicit store register for solution
            this->mRegisterId = EXPLICIT_REGISTER;

         } else if(this->mRegisterId == EXPLICIT_REGISTER)
         {
            // Update explicit term
            for(size_t i = this->mZeroIdx; i < this->mRHSData.size(); i++)
            {
               internal::computeMV(this->mRHSData.at(i), this->mRHSMatrix.at(i), this->mExSolution.at(i));
            }

            // Set ID for solver
            this->mId = TimeSchemeSelector::aIm(this->mStep, this->mStep);

            // Set implicit store register for solution
            this->mRegisterId = IMPLICIT_REGISTER;
         }

         return true;
      }
   }

   template <typename TOperator,typename TData> bool SparseImExRK3RTimestepper<TOperator,TData>::postSolve()
   {
      if(this->mRegisterId == EXPLICIT_REGISTER)
      {
         // Update explicit term
         for(size_t i = this->mZeroIdx; i < this->mRHSData.size(); i++)
         {
            internal::computeSet(this->mExSolution.at(i), this->mSolution.at(i));
         }

         // Loop back to presolve but without new nonlinear term
         this->mHasExplicit = false;

         return true;

      } else if(this->mRegisterId == IMPLICIT_REGISTER)
      {
         // Update implicit term
         for(size_t i = this->mZeroIdx; i < this->mRHSData.size(); i++)
         {
            internal::computeSet(this->mImSolution.at(i), this->mSolution.at(i));
         }

      } else if(this->mRegisterId == SOLUTION_REGISTER)
      {
         // Update intermediate term
         for(size_t i = this->mZeroIdx; i < this->mRHSData.size(); i++)
         {
            internal::computeSet(this->mIntSolution.at(i), this->mSolution.at(i));
         }

         // Loop back to presolve but without new nonlinear term
         this->mHasExplicit = false;

         return true;

      } else if(this->mRegisterId == ERROR_REGISTER)
      {
         // Update error term
         for(size_t i = this->mZeroIdx; i < this->mRHSData.size(); i++)
         {
            internal::computeSet(this->mErrSolution.at(i), this->mSolution.at(i));
         }

         // Loop back to presolve but without new nonlinear term
         this->mHasExplicit = false;

         return true;
      }

      if(this->mStep == 0)
      {
         // Update intermediate solution
         MHDFloat bIm = TimeSchemeSelector::bIm(this->mStep)*this->mDt;
         MHDFloat bEx = TimeSchemeSelector::bEx(this->mStep)*this->mDt;
         for(size_t i = this->mZeroIdx; i < this->mRHSData.size(); i++)
         {
            internal::computeAMXPBYPMZ(this->mIntSolution.at(i), this->mMassMatrix.at(i), bIm, this->mImSolution.at(i), bEx, this->mExSolution.at(i));
         }

         // Embedded lower order scheme solution
         if(TimeSchemeSelector::HAS_EMBEDDED)
         {
            bIm = TimeSchemeSelector::bImErr(this->mStep)*this->mDt;
            bEx = TimeSchemeSelector::bExErr(this->mStep)*this->mDt;
            for(size_t i = this->mZeroIdx; i < this->mRHSData.size(); i++)
            {
               internal::computeAMXPBYPMZ(this->mErrSolution.at(i), this->mMassMatrix.at(i), bIm, this->mImSolution.at(i), bEx, this->mExSolution.at(i));
            }
         }

         // Increase step counter
         this->mStep += 1;

         // Loop back to presolve but without new nonlinear term
         this->mHasExplicit = false;

         return true;

      } else
      {
         // Prepare solution for new nonlinear term
         MHDFloat aIm = TimeSchemeSelector::aIm(this->mStep, this->mStep)*this->mDt;
         for(size_t i = this->mZeroIdx; i < this->mRHSData.size(); i++)
         {
            internal::computeXPAY(this->mSolution.at(i), this->mExSolution.at(i), aIm);
         }

         // Next step will have nonlinear term
         this->mHasExplicit = true;

         return false;
      }
   }

   template <typename TOperator,typename TData> void  SparseImExRK3RTimestepper<TOperator,TData>::updateTimeMatrix(const MHDFloat dt)
   {
      // Update stored timestep
      MHDFloat oldDt = this->mDt;
      this->mDt = dt;

      for(int step = 0; step < TimeSchemeSelector::STEPS; ++step)
      {
         // Update is only required if aIm is not zero
         if(TimeSchemeSelector::aIm(step,step) != 0.0)
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
   }

   template <typename TOperator,typename TData> void SparseImExRK3RTimestepper<TOperator,TData>::initMatrices(const int n)
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

      // Do not reinitialise if work already done by other field
      if(this->mMassMatrix.size() == 0)
      {
         // Reserve space for the RHS matrices
         this->mMassMatrix.reserve(n);

         // Initialise storage for RHS matrices
         for(int i = 0; i < n; ++i)
         {
            // Create storage for LHS matrices
            this->mMassMatrix.push_back(SparseMatrix());
         }
      }
   }

   template <typename TOperator,typename TData> void SparseImExRK3RTimestepper<TOperator,TData>::initSolutions()
   {
      for(size_t i = this->mZeroIdx; i < this->mRHSData.size(); i++)
      {
         internal::computeSet(this->mIntSolution.at(i), this->mSolution.at(i));

         if(TimeSchemeSelector::HAS_EMBEDDED)
         {
            internal::computeSet(this->mErrSolution.at(i), this->mSolution.at(i));
         }
      }
   }

   template <typename TOperator,typename TData> TOperator& SparseImExRK3RTimestepper<TOperator,TData>::rRHSMatrix(const int idx)
   {
      return this->mRHSMatrix.at(idx);
   }

   template <typename TOperator,typename TData> void SparseImExRK3RTimestepper<TOperator,TData>::addStorage(const int rows, const int cols)
   {
      // Assert for non zero rows and columns
      assert(rows > 0);
      assert(cols > 0);

      Solver::SparseLinearSolver<TOperator,TData>::addStorage(rows,cols);

      // Add RK storage
      this->mImSolution.push_back(TData(rows,cols));
      this->mExSolution.push_back(TData(rows,cols));
      this->mTmpSolution.push_back(TData(rows,cols));
      this->mIntSolution.push_back(TData(rows,cols));

      // Initialize storage for embedded scheme
      if(TimeSchemeSelector::HAS_EMBEDDED)
      {
         this->mErrSolution.push_back(TData(rows,cols));
      }
   }

   template <typename TOperator,typename TData> void SparseImExRK3RTimestepper<TOperator,TData>::buildOperators(const int idx, const DecoupledZSparse& opA, const DecoupledZSparse& opB, const DecoupledZSparse& opC, const MHDFloat dt, const int size)
   {
      // Update timestep
      this->mDt = dt;

      // Set explicit matrix
      this->rRHSMatrix(idx).resize(size, size);
      Solver::internal::addOperators(this->rRHSMatrix(idx), 1.0, opA);

      // Set mass matrix
      this->mMassMatrix.at(idx).resize(size, size);
      Solver::internal::addOperators(this->mMassMatrix.at(idx), 1.0, opB);

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

      template <typename TData> inline void computeSet(TData& y, const MHDFloat a, const TData& x)
      {
         y = a*x;
      }

      inline void computeSet(DecoupledZMatrix& y, const MHDFloat a, const DecoupledZMatrix& x)
      {
         y.real() = a*x.real();

         y.imag() = a*x.imag();
      }

      template <typename TData> inline void computeAXPY(TData& y, const MHDFloat a, const TData& x)
      {
         if(a != 0.0)
         {
            y += a*x;
         }
      }

      inline void computeAXPY(DecoupledZMatrix& y, const MHDFloat a, const DecoupledZMatrix& x)
      {
         if(a != 0.0)
         {
            y.real() += a*x.real();

            y.imag() += a*x.imag();
         }
      }

      template <typename TData> inline void computeXPAY(TData& y, const TData& x, const MHDFloat a)
      {
         if(a != 0.0)
         {
            y = x + a*y;

         } else
         {
            computeSet<TData>(y, x);
         }
      }

      inline void computeXPAY(DecoupledZMatrix& y, const DecoupledZMatrix& x, const MHDFloat a)
      {
         if(a != 0.0)
         {
            y.real() = x.real() + a*y.real();

            y.imag() = x.imag() + a*y.imag();

         } else
         {
            computeSet(y, x);
         }
      }

      template <typename TData> inline void computeXPAYPBZ(TData& z, const TData& x, const MHDFloat a, const TData& y, const MHDFloat b)
      {
         if(a == 0.0)
         {
            z = x + b*z;

         } else if(b == 0.0)
         {
            z = x + a*y;

         } else
         {
            z = x + a*y + b*z;
         }
      }

      inline void computeXPAYPBZ(DecoupledZMatrix& z, const DecoupledZMatrix& x, const MHDFloat a, const DecoupledZMatrix& y, const MHDFloat b)
      {
         if(a == 0.0)
         {
            z.real() = x.real() + b*z.real();

            z.imag() = x.imag() + b*z.imag();

         } else if(b == 0.0)
         {
            z.real() = x.real() + a*y.real();

            z.imag() = x.imag() + a*y.imag();

         } else
         {
            z.real() = x.real() + a*y.real() + b*z.real();

            z.imag() = x.imag() + a*y.imag() + b*z.imag();
         }
      }

      template <typename TData> inline void computeAXPBYPZ(TData& z, const MHDFloat a, const TData& x, const MHDFloat b, const TData& y)
      {
         if(a == 0.0)
         {
            z += b*y;

         } else if(b == 0.0)
         {
            z += a*x;

         } else
         {
            z += a*x + b*y;
         }
      }

      inline void computeAXPBYPZ(DecoupledZMatrix& z, const MHDFloat a, const DecoupledZMatrix& x, const MHDFloat b, const DecoupledZMatrix& y)
      {
         if(a == 0.0)
         {
            z.real() += b*y.real();

            z.imag() += b*y.imag();

         } else if(b == 0.0)
         {
            z.real() += a*x.real();

            z.imag() += a*x.imag();

         } else
         {
            z.real() += a*x.real() + b*y.real();

            z.imag() += a*x.imag() + b*y.imag();
         }
      }

      template <typename TData> void computeAMXPY(TData& y, const SparseMatrix& mat, const MHDFloat a, const TData& x)
      {
         if(a != 0.0)
         {
            y += mat*(a*x);
         }
      }

      inline void computeAMXPY(DecoupledZMatrix& y, const SparseMatrix& mat, const MHDFloat a, const DecoupledZMatrix& x)
      {
         if(a != 0.0)
         {
            y.real() += mat*(a*x.real());

            y.imag() += mat*(a*x.imag());
         }
      }

      template <typename TData> void computeAMXPYPBZ(TData& z, const SparseMatrix& mat, const MHDFloat a, const TData& x, const TData& y, const MHDFloat b)
      {
         if(a == 0.0)
         {
            z = y + b*z;

         } else
         {
            z = mat*(a*x) + y + b*z;
         }
      }

      inline void computeAMXPYPBZ(DecoupledZMatrix& z, const SparseMatrix& mat, const MHDFloat a, const DecoupledZMatrix& x, const DecoupledZMatrix& y, const MHDFloat b)
      {
         if(a == 0.0)
         {
            z.real() = y.real() + b*z.real();

            z.imag() = y.imag() + b*z.imag();

         } else
         {
            z.real() = mat*(a*x.real()) + y.real() + b*z.real();

            z.imag() = mat*(a*x.imag()) + y.imag() + b*z.imag();
         }
      }

      template <typename TData> void computeAMXPBYPZ(TData& z, const SparseMatrix& mat, const MHDFloat a, const TData& x, const MHDFloat b, const TData& y)
      {
         if(a == 0.0)
         {
            z += b*y;

         } else if(b == 0.0)
         {
            z += mat*(a*x);

         } else
         {
            z += mat*(a*x) + b*y;
         }
      }

      inline void computeAMXPBYPZ(DecoupledZMatrix& z, const SparseMatrix& mat, const MHDFloat a, const DecoupledZMatrix& x, const MHDFloat b, const DecoupledZMatrix& y)
      {
         if(a == 0.0)
         {
            z.real() += b*y.real();

            z.imag() += b*y.imag();

         } else if(b == 0.0)
         {
            z.real() += mat*(a*x.real());

            z.imag() += mat*(a*x.imag());

         } else
         {
            z.real() += mat*(a*x.real()) + b*y.real();

            z.imag() += mat*(a*x.imag()) + b*y.imag();
         }
      }

      template <typename TData> void computeAMXPBYPMZ(TData& z, const SparseMatrix& mat, const MHDFloat a, const TData& x, const MHDFloat b, const TData& y)
      {
         if(a == 0.0)
         {
            z = b*y + mat*z;

         } else if(b == 0.0)
         {
            z = mat*(a*x + z);

         } else
         {
            z = mat*(a*x + z) + b*y;
         }
      }

      inline void computeAMXPBYPMZ(DecoupledZMatrix& z, const SparseMatrix& mat, const MHDFloat a, const DecoupledZMatrix& x, const MHDFloat b, const DecoupledZMatrix& y)
      {
         if(a == 0.0)
         {
            z.real() = b*y.real() + mat*z.real();

            z.imag() = b*y.imag() + mat*z.imag();

         } else if(b == 0.0)
         {
            z.real() = mat*(a*x.real() + z.real());

            z.imag() = mat*(a*x.imag() + z.imag());

         } else
         {
            z.real() = mat*(a*x.real() + z.real()) + b*y.real();

            z.imag() = mat*(a*x.imag() + z.imag()) + b*y.imag();
         }
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

      template <typename TData> inline void computeError(MHDFloat& err, const TData& x, const TData& y)
      {
         err = std::max(err, ((x.array() - y.array())/(1.0 + x.array().abs())).abs().maxCoeff());
      }

      inline void computeError(MHDFloat& err, const DecoupledZMatrix& x, const DecoupledZMatrix& y)
      {
         err = std::max(err, ((x.real().array() - y.real().array())/(1.0 + x.real().array().abs())).abs().maxCoeff());

         err = std::max(err, ((x.imag().array() - y.imag().array())/(1.0 + x.imag().array().abs())).abs().maxCoeff());
      }

   }
}
}

#endif // SPARSEIMEXRK3RTIMESTEPPER_HPP
