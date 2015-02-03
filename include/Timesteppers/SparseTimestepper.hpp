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
      template <typename TOperator,typename TData> void computeRHSNoMemory(const int step, TData& rRHS, const TOperator& mat, const TData& sol);

      void computeRHSRealNoMemory(const int step, DecoupledZMatrix& rRHS, const SparseMatrix& mat, const DecoupledZMatrix& sol);

      void computeRHSImagNoMemory(const int step, DecoupledZMatrix& rRHS, const SparseMatrix& mat, const DecoupledZMatrix& sol);

      template <typename TOperator,typename TData> void computeRHSNoFieldMemory(const int step, TData& rRHS, const TOperator& mat, const TData& sol, std::vector<TData>& rOldNL, TData& rTmp);

      template <typename TOperator,typename TData> void computeRHSNoNLMemory(const int step, TData& rRHS, const TOperator& mat, const TData& sol, const std::vector<TOperator>& oldMat, std::vector<TData>& rOldSol);

      template <typename TOperator,typename TData> void computeRHSMemory(const int step, TData& rRHS, const TOperator& mat, const TData& sol, const std::vector<TOperator>& oldMat, std::vector<TData>& rOldSol, std::vector<TData>& rOldNL, TData& rTmp);
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
          * @param time    Solver timing with respect to timestepping
          */
         SparseTimestepper(const int start, const SolveTiming::Id time);

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
          * @param lhsCoeff      New coefficient for LHS time dependent part
          * @param rhsCoeff      New coefficient for RHS time dependent part at t_n
          * @param oldRhsCoeff   New coefficient for RHS time dependent part at t_(n-i), i > 0
          * @param step          Timestep scheme substep
          */
         void updateTimeMatrix(const MHDFloat lhsCoeff, const MHDFloat rhsCoeff, const std::vector<MHDFloat>& oldRhsCoeff, const int step);

         /**
          * @brief Set RHS matrix at t_n
          *
          * @param idx Index of the matrix
          */
         TOperator& rRHSMatrix(const int idx);

         /**
          * @brief Set RHS matrix at t_(n-i), i > 0
          *
          * @param i    Time index
          * @param idx  Index of the matrix
          */
         TOperator& rOldRHSMatrix(const int i, const int idx);

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
          * @brief RHS operators of the timestepped equations at t_n
          */
         std::vector<TOperator>   mRHSMatrix;

         /**
          * @brief RHS operators of the timestepped equations at t_(n-i), i > 0
          */
         std::vector<std::vector<TOperator> >   mOldRHSMatrix;

         /**
          * @brief Time dependent part of LHS matrix
          */
         std::vector<TOperator>   mTMatrix;

         /**
          * @brief Storage for old solution at t_(n-i)
          */
         std::vector<std::vector<TData> >  mOldSolution;

         /**
          * @brief Storage for old nonlinear RHS at t_(n-i)
          */
         std::vector<std::vector<TData> >  mOldRHS;

      private:
   };

   template <typename TOperator,typename TData> SparseTimestepper<TOperator,TData>::SparseTimestepper(const int start, const SolveTiming::Id time)
      : Solver::SparseLinearSolver<TOperator,TData>(start, time)
   {
   }

   template <typename TOperator,typename TData> SparseTimestepper<TOperator,TData>::~SparseTimestepper()
   {
   }

   template <typename TOperator,typename TData> void SparseTimestepper<TOperator,TData>::computeRHS(const int step)
   {
      int start = step*this->nSystem();

      if(IntegratorSelector::FIELD_MEMORY == 0 && IntegratorSelector::NONLINEAR_MEMORY == 0)
      {
         for(size_t i = this->mZeroIdx; i < this->mRHSData.size(); i++)
         {
            internal::computeRHSNoMemory<TOperator,TData>(step, this->mRHSData.at(i), this->mRHSMatrix.at(i+start), this->mSolution.at(i));
         }

      } else if(IntegratorSelector::FIELD_MEMORY == 0)
      {
         TData tmp;

         for(size_t i = this->mZeroIdx; i < this->mRHSData.size(); i++)
         {
            internal::computeRHSNoFieldMemory<TOperator,TData>(step, this->mRHSData.at(i), this->mRHSMatrix.at(i+start), this->mSolution.at(i), this->mOldRHS.at(i), tmp);
         }

      } else if(IntegratorSelector::NONLINEAR_MEMORY == 0)
      {
         for(size_t i = this->mZeroIdx; i < this->mRHSData.size(); i++)
         {
            internal::computeRHSNoNLMemory<TOperator,TData>(step, this->mRHSData.at(i), this->mRHSMatrix.at(i+start), this->mSolution.at(i), this->mOldRHSMatrix.at(i+start), this->mOldSolution.at(i));
         }

      } else
      {
         TData tmp;

         for(size_t i = this->mZeroIdx; i < this->mRHSData.size(); i++)
         {
            internal::computeRHSMemory<TOperator,TData>(step, this->mRHSData.at(i), this->mRHSMatrix.at(i+start), this->mSolution.at(i), this->mOldRHSMatrix.at(i+start), this->mOldSolution.at(i), this->mOldRHS.at(i), tmp);
         }
      }
   }

   template <typename TOperator,typename TData> void  SparseTimestepper<TOperator,TData>::updateTimeMatrix(const MHDFloat lhsCoeff, const MHDFloat rhsCoeff, const std::vector<MHDFloat>& oldRhsCoeff, const int step)
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
         std::vector<size_t> oldRhsJ;
         for(int r = 0; r < IntegratorSelector::FIELD_MEMORY; r++)
         {
            oldRhsJ.push_back(0);
         }
         for (size_t k=0; k< static_cast<size_t>(this->mTMatrix.at(i).outerSize()); ++k)
         {
            typename TOperator::InnerIterator lhsIt(this->mLHSMatrix.at(start+i),lhsJ);
            typename TOperator::InnerIterator rhsIt(this->mRHSMatrix.at(start+i),rhsJ);
            std::vector<SharedPtrMacro<typename TOperator::InnerIterator> > oldRhsIt;
            for(int r = 0; r < IntegratorSelector::FIELD_MEMORY; r++)
            {
               SharedPtrMacro<typename TOperator::InnerIterator> spTmpIt(new typename TOperator::InnerIterator(this->mOldRHSMatrix.at(r).at(start+i),oldRhsJ.at(r)));
               oldRhsIt.push_back(spTmpIt);
            }
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

                  // Update RHS matrix at t_(n-i), i > 0
                  for(int r = 0; r < IntegratorSelector::FIELD_MEMORY; r++)
                  {
                     if(timeIt.col() == oldRhsIt.at(r)->col())
                     {
                        if(timeIt.row() == oldRhsIt.at(r)->row())
                        {
                           // Update values
                           oldRhsIt.at(r)->valueRef() += oldRhsCoeff.at(r)*timeIt.value();

                           // Update RHS iterators and counters
                           ++(*oldRhsIt.at(r));
                           if(!(*oldRhsIt.at(r)))
                           {
                              oldRhsJ.at(r)++;
                           }
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
         // Reserve space for the RHS matrices at t_n
         this->mRHSMatrix.reserve(n);

         // Initialise storage at t_n
         for(int i = 0; i < n; ++i)
         {
            // Create storage for LHS matrices
            this->mRHSMatrix.push_back(TOperator());
         }

         // Reserve space for the RHS matrices at t_(n-i), i > 0 
         if(IntegratorSelector::FIELD_MEMORY > 0)
         {
            this->mOldRHSMatrix.reserve(n);

            // Initialise storage at t_(n-i), i > 0
            for(int i = 0; i < n; ++i)
            {
               this->mOldRHSMatrix.push_back(std::vector<TOperator>());
               this->mOldRHSMatrix.back().reserve(IntegratorSelector::FIELD_MEMORY);

               // Handle RHS matrices at t_(n-i), i > 0
               for(int j = 0; j < IntegratorSelector::FIELD_MEMORY; j++)
               {
                  // Create storage for LHS matrices
                  this->mOldRHSMatrix.back().push_back(TOperator());
               }
            }

         }
      }
   }

   template <typename TOperator,typename TData> TOperator& SparseTimestepper<TOperator,TData>::rRHSMatrix(const int idx)
   {
      return this->mRHSMatrix.at(idx);
   }

   template <typename TOperator,typename TData> TOperator& SparseTimestepper<TOperator,TData>::rOldRHSMatrix(const int i, const int idx)
   {
      return this->mOldRHSMatrix.at(idx).at(i);
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

      // Add storage for RHS at t_(n-i), i > 0
      if(IntegratorSelector::NONLINEAR_MEMORY > 0)
      {
         this->mOldRHS.push_back(std::vector<TData>());
         for(int j = 0; j < IntegratorSelector::NONLINEAR_MEMORY; j++)
         {
            this->mOldRHS.back().push_back(TData(rows,cols));
            this->mOldRHS.back().back().setZero();
         }
      }

      // Add storage for solution at t_(n-i), i > 0
      if(IntegratorSelector::FIELD_MEMORY > 0)
      {
         this->mOldSolution.push_back(std::vector<TData>());
         for(int j = 0; j < IntegratorSelector::FIELD_MEMORY; j++)
         {
            this->mOldSolution.back().push_back(TData(rows,cols));
            this->mOldSolution.back().back().setZero();
         }
      }

      // Add storage for the time matrix
      this->mTMatrix.push_back(TOperator());
   }

   namespace internal
   {
      template <typename TOperator,typename TData> inline void computeRHSNoMemory(const int step, TData& rRHS, const TOperator& mat, const TData& sol)
      {
         rRHS = mat*sol + IntegratorSelector::rhsN(0, step)*rRHS;
      }

      template <typename TOperator,typename TData> inline void computeRHSNoFieldMemory(const int step, TData& rRHS, const TOperator& mat, const TData& sol, std::vector<TData>& rOldNL, TData& rTmp)
      {
         // Store RHS at t_n
         rTmp = rRHS;

         // Compute RHS part from solution and nonlinear term at t_n
         computeRHSNoMemory(step, rRHS, mat, sol);

         // Compute RHS part from nonlinear terms at t_(n-i), i > 0
         for(int n = 0; n < IntegratorSelector::nonlinearMemory(step); n++)
         {
            rRHS += IntegratorSelector::rhsN(n+1, step)*rOldNL.at(n);
         }

         // Shift old RHS by one
         for(int n = IntegratorSelector::NONLINEAR_MEMORY-1; n > 0; n--)
         {
            rOldNL.at(n) = rOldNL.at(n-1); 
         }
         rOldNL.at(0) = rTmp;
      }

      template <typename TOperator,typename TData> inline void computeRHSNoNLMemory(const int step, TData& rRHS, const TOperator& mat, const TData& sol, const std::vector<TOperator>& oldMat, std::vector<TData>& rOldSol)
      {
         // Compute RHS part from solution and nonlinear term at t_n
         computeRHSNoMemory(step, rRHS, mat, sol);

         // Compute RHS part from solution at t_(n-i), i > 0
         for(int i = 0; i < IntegratorSelector::fieldMemory(step); i++)
         {
            rRHS += oldMat.at(i)*rOldSol.at(i);
         }

         // Shift old solution storage by one
         for(int n = IntegratorSelector::FIELD_MEMORY-1; n > 0; n--)
         {
            rOldSol.at(n) = rOldSol.at(n-1); 
         }
         rOldSol.at(0) = sol;
      }

      template <typename TOperator,typename TData> inline void computeRHSMemory(const int step, TData& rRHS, const TOperator& mat, const TData& sol, const std::vector<TOperator>& oldMat, std::vector<TData>& rOldSol, std::vector<TData>& rOldNL, TData& rTmp)
      {
         // Store RHS at t_n
         rTmp = rRHS;

         // Compute RHS part from solution and nonlinear term at t_n
         computeRHSNoMemory(step, rRHS, mat, sol);

         // Compute RHS part from solution at t_(n-i), i > 0
         for(int i = 0; i < IntegratorSelector::fieldMemory(step); i++)
         {
            rRHS += oldMat.at(i)*rOldSol.at(i);
         }

         // Shift old solution storage by one
         for(int n = IntegratorSelector::FIELD_MEMORY-1; n > 0; n--)
         {
            rOldSol.at(n) = rOldSol.at(n-1); 
         }
         rOldSol.at(0) = sol;

         // Compute RHS part from nonlinear terms at t_(n-i), i > 0
         for(int n = 0; n < IntegratorSelector::nonlinearMemory(step); n++)
         {
            rRHS += IntegratorSelector::rhsN(n+1, step)*rOldNL.at(n);
         }

         // Shift old storage by one
         for(int n = IntegratorSelector::NONLINEAR_MEMORY-1; n > 0; n--)
         {
            rOldNL.at(n) = rOldNL.at(n-1); 
         }
         rOldNL.at(0) = rTmp;
      }

      inline void computeRHSRealNoMemory(const int step, DecoupledZMatrix& rRHS, const SparseMatrix& mat, const DecoupledZMatrix& sol)
      {
         rRHS.real() = mat*sol.real() + IntegratorSelector::rhsN(0, step)*rRHS.real();
      }

      inline void computeRHSImagNoMemory(const int step, DecoupledZMatrix& rRHS, const SparseMatrix& mat, const DecoupledZMatrix& sol)
      {
         rRHS.imag() = mat*sol.imag() + IntegratorSelector::rhsN(0, step)*rRHS.imag();
      }

      template <> inline void computeRHSNoMemory<SparseMatrix,DecoupledZMatrix>(const int step, DecoupledZMatrix& rRHS, const SparseMatrix& mat, const DecoupledZMatrix& sol)
      {
         computeRHSRealNoMemory(step, rRHS, mat, sol);

         computeRHSImagNoMemory(step, rRHS, mat, sol);
      }

      template <> inline void computeRHSNoFieldMemory<SparseMatrix,DecoupledZMatrix>(const int step, DecoupledZMatrix& rRHS, const SparseMatrix& mat, const DecoupledZMatrix& sol, std::vector<DecoupledZMatrix>& rOldNL, DecoupledZMatrix& rTmp)
      {
         //
         // rTmp is only a temporary variable, using only real reduces memory usage
         //
         
         // Real part
         
         // Store RHS at t_n
         rTmp.real() = rRHS.real();

         // Compute RHS part from solution and nonlinear term at t_n
         computeRHSRealNoMemory(step, rRHS, mat, sol);

         // Compute RHS part from nonlinear terms at t_(n-i), i > 0
         for(int n = 0; n < IntegratorSelector::nonlinearMemory(step); n++)
         {
            rRHS.real() += IntegratorSelector::rhsN(n+1, step)*rOldNL.at(n).real();
         }

         // Shift old storage by one
         for(int n = IntegratorSelector::NONLINEAR_MEMORY-1; n > 0; n--)
         {
            rOldNL.at(n).real() = rOldNL.at(n-1).real(); 
         }
         rOldNL.at(0).real() = rTmp.real();

         // Imaginary part
         
         // Store RHS at t_n
         rTmp.real() = rRHS.imag();

         // Compute RHS part from solution and nonlinear term at t_n
         computeRHSImagNoMemory(step, rRHS, mat, sol);

         // Compute RHS part from nonlinear terms at t_(n-i), i > 0
         for(int n = 0; n < IntegratorSelector::nonlinearMemory(step); n++)
         {
            rRHS.imag() += IntegratorSelector::rhsN(n+1, step)*rOldNL.at(n).imag();
         }

         // Shift old storage by one
         for(int n = IntegratorSelector::NONLINEAR_MEMORY-1; n > 0; n--)
         {
            rOldNL.at(n).imag() = rOldNL.at(n-1).imag(); 
         }
         rOldNL.at(0).imag() = rTmp.real();
      }

      template <> inline void computeRHSNoNLMemory<SparseMatrix,DecoupledZMatrix>(const int step, DecoupledZMatrix& rRHS, const SparseMatrix& mat, const DecoupledZMatrix& sol, const std::vector<SparseMatrix>& oldMat, std::vector<DecoupledZMatrix>& rOldSol)
      {
         // Compute RHS part from solution and nonlinear term at t_n
         computeRHSNoMemory(step, rRHS, mat, sol);

         // Real part

         // Compute RHS part from solution at t_(n-i), i > 0
         for(int i = 0; i < IntegratorSelector::fieldMemory(step); i++)
         {
            rRHS.real() += oldMat.at(i)*rOldSol.at(i).real();
         }

         // Shift old solution storage by one
         for(int n = IntegratorSelector::FIELD_MEMORY-1; n > 0; n--)
         {
            rOldSol.at(n).real() = rOldSol.at(n-1).real(); 
         }
         rOldSol.at(0).real() = sol.real();

         // Imaginary part

         // Compute RHS part from solution at t_(n-i), i > 0
         for(int i = 0; i < IntegratorSelector::fieldMemory(step); i++)
         {
            rRHS.imag() += oldMat.at(i)*rOldSol.at(i).imag();
         }

         // Shift old solution storage by one
         for(int n = IntegratorSelector::FIELD_MEMORY-1; n > 0; n--)
         {
            rOldSol.at(n).imag() = rOldSol.at(n-1).imag(); 
         }
         rOldSol.at(0).imag() = sol.imag();
      }

      template <> inline void computeRHSMemory<SparseMatrix,DecoupledZMatrix>(const int step, DecoupledZMatrix& rRHS, const SparseMatrix& mat, const DecoupledZMatrix& sol, const std::vector<SparseMatrix>& oldMat, std::vector<DecoupledZMatrix>& rOldSol, std::vector<DecoupledZMatrix>& rOldNL, DecoupledZMatrix& rTmp)
      {
         //
         // rTmp is only a temporary variable, using only real reduces memory usage
         //
         
         // Real part
         
         // Store RHS at t_n
         rTmp.real() = rRHS.real();

         // Compute RHS part from solution and nonlinear term at t_n
         computeRHSRealNoMemory(step, rRHS, mat, sol);

         // Compute RHS part from solution at t_(n-i), i > 0
         for(int i = 0; i < IntegratorSelector::fieldMemory(step); i++)
         {
            rRHS.real() += oldMat.at(i)*rOldSol.at(i).real();
         }

         // Shift old solution storage by one
         for(int n = IntegratorSelector::FIELD_MEMORY-1; n > 0; n--)
         {
            rOldSol.at(n).real() = rOldSol.at(n-1).real(); 
         }
         rOldSol.at(0).real() = sol.real();

         // Compute RHS part from nonlinear terms at t_(n-i), i > 0
         for(int n = 0; n < IntegratorSelector::nonlinearMemory(step); n++)
         {
            rRHS.real() += IntegratorSelector::rhsN(n+1, step)*rOldNL.at(n).real();
         }

         // Shift old storage by one
         for(int n = IntegratorSelector::NONLINEAR_MEMORY-1; n > 0; n--)
         {
            rOldNL.at(n).real() = rOldNL.at(n-1).real(); 
         }
         rOldNL.at(0).real() = rTmp.real();

         // Imaginary part
         
         // Store RHS at t_n
         rTmp.real() = rRHS.imag();

         // Compute RHS part from solution and nonlinear term at t_n
         computeRHSImagNoMemory(step, rRHS, mat, sol);

         // Compute RHS part from solution at t_(n-i), i > 0
         for(int i = 0; i < IntegratorSelector::fieldMemory(step); i++)
         {
            rRHS.imag() += oldMat.at(i)*rOldSol.at(i).imag();
         }

         // Shift old solution storage by one
         for(int n = IntegratorSelector::FIELD_MEMORY-1; n > 0; n--)
         {
            rOldSol.at(n).imag() = rOldSol.at(n-1).imag(); 
         }
         rOldSol.at(0).imag() = sol.imag();

         // Compute RHS part from nonlinear terms at t_(n-i), i > 0
         for(int n = 0; n < IntegratorSelector::nonlinearMemory(step); n++)
         {
            rRHS.imag() += IntegratorSelector::rhsN(n+1, step)*rOldNL.at(n).imag();
         }

         // Shift old storage by one
         for(int n = IntegratorSelector::NONLINEAR_MEMORY-1; n > 0; n--)
         {
            rOldNL.at(n).imag() = rOldNL.at(n-1).imag(); 
         }
         rOldNL.at(0).imag() = rTmp.real();
      }
   }
}
}

#endif // SPARSETIMESTEPPER_HPP
