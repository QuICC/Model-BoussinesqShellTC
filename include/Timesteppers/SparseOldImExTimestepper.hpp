/** 
 * @file SparseOldImExTimestepper.hpp
 * @brief Implementation of a templated (coupled) equation timestepper
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef SPARSEOLDIMEXTIMESTEPPER_HPP
#define SPARSEOLDIMEXTIMESTEPPER_HPP

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

namespace QuICC {

namespace Timestep {

   namespace internal
   {
      template <typename TOperator,typename TData> void computeMVPAY(TData& y, const TOperator& A, const TData& x, const MHDFloat a);

      void computeMVPAY(DecoupledZMatrix& y, const SparseMatrix& A, const DecoupledZMatrix& x, const MHDFloat a);


      template <typename TOperator,typename TData> void computeRHSNoFieldMemory(const int step, TData& rRHS, const TOperator& mat, const TData& sol, std::vector<TData>& rOldNL, TData& rTmp);

      template <typename TOperator,typename TData> void computeRHSNoNLMemory(const int step, TData& rRHS, const TOperator& mat, const TData& sol, const std::vector<TOperator>& oldMat, std::vector<TData>& rOldSol);

      template <typename TOperator,typename TData> void computeRHSMemory(const int step, TData& rRHS, const TOperator& mat, const TData& sol, const std::vector<TOperator>& oldMat, std::vector<TData>& rOldSol, std::vector<TData>& rOldNL, TData& rTmp);
   }

   /**
    * @brief Implementation of a templated (coupled) equation timestepper
    */
   template <typename TOperator,typename TData> class SparseOldImExTimestepper: public Solver::SparseLinearSolver<TOperator,TData>
   {
      public:
         /**
          * @brief Constructor
          *
          * @param start   Starting index (for example without m=0)
          * @param time    Solver timing with respect to timestepping
          */
         SparseOldImExTimestepper(const int start, const SolveTiming::Id time);

         /**
          * @brief Destructor
          */
         virtual ~SparseOldImExTimestepper();

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
          * @brief Compute the RHS of the linear systems
          *
          * @param step    Current substep
          */
         void computeRHS(const int step);

         /**
          * @brief Prepare fields for implicit solve
          */
         virtual bool preSolve();

         /**
          * @brief Update the LHS matrix with new timedependence
          */
         void updateTimeMatrix(const MHDFloat dt);

         /**
          * @brief Set LHS matrix
          *
          * @param idx Index of the matrix
          */
         TOperator& rRHSMatrix(const MHDFloat id, const int idx);

         /**
          * @brief Set RHS matrix at t_(n-i), i > 0
          *
          * @param i    Time index
          * @param idx  Index of the matrix
          */
         std::vector<TOperator>& rOldRHSMatrices(const MHDFloat id, const int idx);

         /**
          * @brief Set RHS matrix at t_(n-i), i > 0
          *
          * @param i    Time index
          * @param idx  Index of the matrix
          */
         TOperator& rOldRHSMatrix(const MHDFloat id, const int idx, const int i);

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

         /**
          * @brief Build the scheme operators
          *
          * @param idx  Solver index
          * @param opA  A operator
          * @param opB  B operator
          * @param opC  Boundary operator
          * @param opD  Inhomogeneous boundary operator
          */
         void buildOperators(const int idx, const DecoupledZSparse& opA, const DecoupledZSparse& opB, const DecoupledZSparse& opC, const DecoupledZSparse& opD, const MHDFloat dt, const int size);

         /**
          * @brief Error of computation
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
          * @brief Initialise the solver matrices storage
          *
          * @param n Size of matrices
          */
         void initMatrices(const MHDFloat lhsId, const MHDFloat rhsId, const int n);

         /**
          * @brief Current substep
          */
         int mStep;

         /**
          * @brief Current timestep
          */
         MHDFloat mDt;

         /**
          * @brief RHS operators of the timestepped equations at t_n
          */
         std::map<MHDFloat, std::vector<TOperator> >  mRHSMatrix;

         /**
          * @brief RHS operators of the timestepped equations at t_(n-i), i > 0
          */
         std::map<MHDFloat, std::vector<std::vector<TOperator> > >   mOldRHSMatrix;

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

   template <typename TOperator,typename TData> SparseOldImExTimestepper<TOperator,TData>::SparseOldImExTimestepper(const int start, const SolveTiming::Id time)
      : Solver::SparseLinearSolver<TOperator,TData>(start, time), mStep(0), mDt(-1.0)
   {
   }

   template <typename TOperator,typename TData> SparseOldImExTimestepper<TOperator,TData>::~SparseOldImExTimestepper()
   {
   }

   template <typename TOperator,typename TData> MHDFloat SparseOldImExTimestepper<TOperator,TData>::error() const
   {
      return -1.0;
   }

   template <typename TOperator,typename TData> bool SparseOldImExTimestepper<TOperator,TData>::finished()
   {
      this->mStep = (this->mStep + 1) % (TimeSchemeSelector::STEPS);

      return (this->mStep == 0);
   }

   template <typename TOperator,typename TData> MHDFloat SparseOldImExTimestepper<TOperator,TData>::stepFraction() const
   {
      return TimeSchemeSelector::cEx(this->mStep);
   }

   template <typename TOperator,typename TData> bool SparseOldImExTimestepper<TOperator,TData>::preSolve()
   {
      this->mId = TimeSchemeSelector::lhsT(this->mStep);
      MHDFloat rhsId = TimeSchemeSelector::rhsT(0,this->mStep);

      if(TimeSchemeSelector::FIELD_MEMORY == 0 && TimeSchemeSelector::NONLINEAR_MEMORY == 0)
      {
         for(size_t i = this->mZeroIdx; i < this->mRHSData.size(); i++)
         {
            internal::computeMVPAY<TOperator,TData>(this->mRHSData.at(i), this->rRHSMatrix(rhsId,i), this->mSolution.at(i), TimeSchemeSelector::rhsN(0, this->mStep));
         }

      } else if(TimeSchemeSelector::FIELD_MEMORY == 0)
      {
         TData tmp;

         for(size_t i = this->mZeroIdx; i < this->mRHSData.size(); i++)
         {
            internal::computeRHSNoFieldMemory<TOperator,TData>(this->mStep, this->mRHSData.at(i), this->rRHSMatrix(rhsId,i), this->mSolution.at(i), this->mOldRHS.at(i), tmp);
         }

      } else if(TimeSchemeSelector::NONLINEAR_MEMORY == 0)
      {
         for(size_t i = this->mZeroIdx; i < this->mRHSData.size(); i++)
         {
            internal::computeRHSNoNLMemory<TOperator,TData>(this->mStep, this->mRHSData.at(i), this->rRHSMatrix(rhsId,i), this->mSolution.at(i), this->rOldRHSMatrices(rhsId,i), this->mOldSolution.at(i));
         }

      } else
      {
         TData tmp;

         for(size_t i = this->mZeroIdx; i < this->mRHSData.size(); i++)
         {
            internal::computeRHSMemory<TOperator,TData>(this->mStep, this->mRHSData.at(i), this->rRHSMatrix(rhsId,i), this->mSolution.at(i), this->rOldRHSMatrices(rhsId,i), this->mOldSolution.at(i), this->mOldRHS.at(i), tmp);
         }
      }

      // Include inhomogeneous boundary conditions
      this->addInhomogeneous();

      return true;
   }

   template <typename TOperator,typename TData> void  SparseOldImExTimestepper<TOperator,TData>::updateTimeMatrix(const MHDFloat dt)
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
            size_t nnz = this->mTMatrix.at(i).nonZeros();

            // Update LHS and RHS matrices
            size_t lhsJ = 0;
            size_t rhsJ = 0;
            std::vector<size_t> oldRhsJ;
            for(int r = 0; r < TimeSchemeSelector::FIELD_MEMORY; r++)
            {
               oldRhsJ.push_back(0);
            }
            for (size_t k=0; k< static_cast<size_t>(this->mTMatrix.at(i).outerSize()); ++k)
            {
               typename TOperator::InnerIterator lhsIt(this->rLHSMatrix(TimeSchemeSelector::lhsT(step), i),lhsJ);
               typename TOperator::InnerIterator rhsIt(this->rRHSMatrix(TimeSchemeSelector::rhsT(0,step), i),rhsJ);
               std::vector<SharedPtrMacro<typename TOperator::InnerIterator> > oldRhsIt;
               for(int r = 0; r < TimeSchemeSelector::FIELD_MEMORY; r++)
               {
                  SharedPtrMacro<typename TOperator::InnerIterator> spTmpIt(new typename TOperator::InnerIterator(this->rOldRHSMatrix(TimeSchemeSelector::rhsT(0,step),i,r),oldRhsJ.at(r)));
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
                           lhsIt.valueRef() += TimeSchemeSelector::lhsT(step)*(1.0/oldDt - 1.0/this->mDt)*timeIt.value();

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
                           rhsIt.valueRef() += TimeSchemeSelector::rhsT(0,step)*(1.0/oldDt - 1.0/this->mDt)*timeIt.value();

                           // Update RHS iterators and counters
                           ++rhsIt;
                           if(!rhsIt)
                           {
                              rhsJ++;
                           }
                        }
                     }

                     // Update RHS matrix at t_(n-i), i > 0
                     for(int r = 0; r < TimeSchemeSelector::FIELD_MEMORY; r++)
                     {
                        if(timeIt.col() == oldRhsIt.at(r)->col())
                        {
                           if(timeIt.row() == oldRhsIt.at(r)->row())
                           {
                              // Update values
                              oldRhsIt.at(r)->valueRef() += TimeSchemeSelector::rhsT(r+1,step)*(1.0/oldDt - 1.0/this->mDt)*timeIt.value();

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
   }

   template <typename TOperator,typename TData> void SparseOldImExTimestepper<TOperator,TData>::initMatrices(const int n)
   {
      // Initialise base matrices
      for(int i = 0; i < TimeSchemeSelector::STEPS; i++)
      {
         this->initMatrices(TimeSchemeSelector::lhsT(i), TimeSchemeSelector::rhsT(0,i), n);
      }
   }

   template <typename TOperator,typename TData> void SparseOldImExTimestepper<TOperator,TData>::initSolutions()
   {
      for(size_t k = 0; k < this->mOldSolution.size(); ++k)
      {
         for(size_t j = 0; j < this->mOldSolution.at(k).size(); ++j)
         {
            this->mOldSolution.at(k).at(j) = this->mSolution.at(k);
         }
      }
   }

   template <typename TOperator,typename TData> TOperator& SparseOldImExTimestepper<TOperator,TData>::rRHSMatrix(const MHDFloat id, const int idx)
   {
      return this->mRHSMatrix.find(id)->second.at(idx);
   }

   template <typename TOperator,typename TData> std::vector<TOperator>& SparseOldImExTimestepper<TOperator,TData>::rOldRHSMatrices(const MHDFloat id, const int idx)
   {
      return this->mOldRHSMatrix.find(id)->second.at(idx);
   }

   template <typename TOperator,typename TData> TOperator& SparseOldImExTimestepper<TOperator,TData>::rOldRHSMatrix(const MHDFloat id, const int idx, const int i)
   {
      return this->mOldRHSMatrix.find(id)->second.at(idx).at(i);
   }

   template <typename TOperator,typename TData> TOperator& SparseOldImExTimestepper<TOperator,TData>::rTMatrix(const int idx)
   {
      return this->mTMatrix.at(idx);
   }

   template <typename TOperator,typename TData> void SparseOldImExTimestepper<TOperator,TData>::addStorage(const int rows, const int cols)
   {
      // Assert for non zero rows and columns
      assert(rows > 0);
      assert(cols > 0);

      Solver::SparseLinearSolver<TOperator,TData>::addStorage(rows,cols);

      // Add storage for RHS at t_(n-i), i > 0
      if(TimeSchemeSelector::NONLINEAR_MEMORY > 0)
      {
         this->mOldRHS.push_back(std::vector<TData>());
         for(int j = 0; j < TimeSchemeSelector::NONLINEAR_MEMORY; j++)
         {
            this->mOldRHS.back().push_back(TData(rows,cols));
            this->mOldRHS.back().back().setZero();
         }
      }

      // Add storage for solution at t_(n-i), i > 0
      if(TimeSchemeSelector::FIELD_MEMORY > 0)
      {
         this->mOldSolution.push_back(std::vector<TData>());
         for(int j = 0; j < TimeSchemeSelector::FIELD_MEMORY; j++)
         {
            this->mOldSolution.back().push_back(TData(rows,cols));
            this->mOldSolution.back().back().setZero();
         }
      }

      // Add storage for the time matrix
      this->mTMatrix.push_back(TOperator());
   }

   template <typename TOperator,typename TData> void SparseOldImExTimestepper<TOperator,TData>::buildOperators(const int idx, const DecoupledZSparse& opA, const DecoupledZSparse& opB, const DecoupledZSparse& opC, const DecoupledZSparse& opD, const MHDFloat dt, const int size)
   {
      // Update timestep
      this->mDt = dt;

      // Set time matrix
      this->rTMatrix(idx).resize(size, size);
      Solver::internal::addOperators(this->rTMatrix(idx), 1.0, opB);

      // Set implicit and explicit matrices
      for(int i = 0; i < TimeSchemeSelector::STEPS; ++i)
      {
         // Set implicit matrix
         this->rLHSMatrix(TimeSchemeSelector::lhsT(i), idx).resize(size, size);
         Solver::internal::addOperators(this->rLHSMatrix(TimeSchemeSelector::lhsT(i), idx), TimeSchemeSelector::lhsL(i), opA);
         Solver::internal::addOperators(this->rLHSMatrix(TimeSchemeSelector::lhsT(i), idx), -TimeSchemeSelector::lhsT(i)/this->mDt, opB);
         Solver::internal::addOperators(this->rLHSMatrix(TimeSchemeSelector::lhsT(i), idx), 1.0, opC);

         // Set explicit matrix
         this->rRHSMatrix(TimeSchemeSelector::rhsT(0,i), idx).resize(size, size);
         this->rRHSMatrix(TimeSchemeSelector::rhsT(0,i), idx).resize(size, size);
         Solver::internal::addOperators(this->rRHSMatrix(TimeSchemeSelector::rhsT(0,i), idx), -TimeSchemeSelector::rhsL(0,i), opA);
         Solver::internal::addOperators(this->rRHSMatrix(TimeSchemeSelector::rhsT(0,i), idx), -TimeSchemeSelector::rhsT(0,i)/this->mDt, opB);

         // Set explicit matrix for t_(n-i)
         for(int j = 0; j < TimeSchemeSelector::FIELD_MEMORY; j++)
         {
            this->rOldRHSMatrix(TimeSchemeSelector::rhsT(0,i), idx, j).resize(size, size);
            this->rOldRHSMatrix(TimeSchemeSelector::rhsT(0,i), idx, j).resize(size, size);
            Solver::internal::addOperators(this->rOldRHSMatrix(TimeSchemeSelector::rhsT(0,i), idx, j), -TimeSchemeSelector::rhsL(j+1,i), opA);
            Solver::internal::addOperators(this->rOldRHSMatrix(TimeSchemeSelector::rhsT(0,i), idx, j), -TimeSchemeSelector::rhsT(j+1,i)/this->mDt, opB);
         }
      }

      // Set inhomogeneous boundary value
      this->setInhomogeneous(idx, opD);
   }

   template <typename TOperator,typename TData> void SparseOldImExTimestepper<TOperator,TData>::initMatrices(const MHDFloat lhsId, const MHDFloat rhsId, const int n)
   {
      Solver::SparseLinearSolver<TOperator,TData>::initMatrices(lhsId, n);

      // Do not reinitialise if work already done by other field
      if(this->mRHSMatrix.count(rhsId) == 0 || this->mRHSMatrix.find(rhsId)->second.size() == 0)
      {  
         if(this->mRHSMatrix.count(rhsId) == 0)
         {
            this->mRHSMatrix.insert(std::make_pair(rhsId, std::vector<TOperator>()));
         }

         typename std::map<MHDFloat,std::vector<TOperator> >::iterator it = this->mRHSMatrix.find(rhsId);

         // Reserve space for the RHS matrices
         it->second.reserve(n);

         // Initialise storage
         for(int i = 0; i < n; ++i)
         {
            // Create storage for RHS matrices
            it->second.push_back(TOperator());
         }
      }

      // Reserve space for the RHS matrices at t_(n-i), i > 0 
      if((this->mOldRHSMatrix.count(rhsId) == 0 || this->mOldRHSMatrix.find(rhsId)->second.size() == 0) && TimeSchemeSelector::FIELD_MEMORY > 0)
      {
         if(this->mOldRHSMatrix.count(rhsId) == 0)
         {
            this->mOldRHSMatrix.insert(std::make_pair(rhsId, std::vector<std::vector<TOperator> >()));
         }

         typename std::map<MHDFloat,std::vector<std::vector<TOperator> > >::iterator it = this->mOldRHSMatrix.find(rhsId);

         // Reserve space for the RHS matrices
         it->second.reserve(n);

         // Initialise storage
         for(int i = 0; i < n; ++i)
         {
            // Create storage for RHS matrices
            it->second.push_back(std::vector<TOperator>());
            it->second.back().reserve(TimeSchemeSelector::FIELD_MEMORY);

            // Handle RHS matrices at t_(n-i), i > 0
            for(int j = 0; j < TimeSchemeSelector::FIELD_MEMORY; j++)
            {
               // Create storage for LHS matrices
               it->second.back().push_back(TOperator());
            }
         }
      }
   }

   namespace internal
   {
      template <typename TOperator,typename TData> inline void computeMVPAY(TData& y, const TOperator& A, const TData& x, const MHDFloat a)
      {
         y = A*x + a*y;
      }

      template <typename TOperator,typename TData> inline void computeRHSNoFieldMemory(const int step, TData& rRHS, const TOperator& mat, const TData& sol, std::vector<TData>& rOldNL, TData& rTmp)
      {
         // Store RHS at t_n
         rTmp = rRHS;

         // Compute RHS part from solution and nonlinear term at t_n
         computeMVPAY<TOperator,TData>(rRHS, mat, sol, TimeSchemeSelector::rhsN(0, step));

         // Compute RHS part from nonlinear terms at t_(n-i), i > 0
         for(int n = 0; n < TimeSchemeSelector::nonlinearMemory(step); n++)
         {
            rRHS += TimeSchemeSelector::rhsN(n+1, step)*rOldNL.at(n);
         }

         // Shift old RHS by one
         for(int n = TimeSchemeSelector::NONLINEAR_MEMORY-1; n > 0; n--)
         {
            rOldNL.at(n) = rOldNL.at(n-1); 
         }
         rOldNL.at(0) = rTmp;
      }

      template <typename TOperator,typename TData> inline void computeRHSNoNLMemory(const int step, TData& rRHS, const TOperator& mat, const TData& sol, const std::vector<TOperator>& oldMat, std::vector<TData>& rOldSol)
      {
         // Compute RHS part from solution and nonlinear term at t_n
         computeMVPAY<TOperator,TData>(rRHS, mat, sol, TimeSchemeSelector::rhsN(0, step));

         // Compute RHS part from solution at t_(n-i), i > 0
         for(int i = 0; i < TimeSchemeSelector::fieldMemory(step); i++)
         {
            rRHS += oldMat.at(i)*rOldSol.at(i);
         }

         // Shift old solution storage by one
         for(int n = TimeSchemeSelector::FIELD_MEMORY-1; n > 0; n--)
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
         computeMVPAY<TOperator,TData>(rRHS, mat, sol, TimeSchemeSelector::rhsN(0, step));

         // Compute RHS part from solution at t_(n-i), i > 0
         for(int i = 0; i < TimeSchemeSelector::fieldMemory(step); i++)
         {
            rRHS += oldMat.at(i)*rOldSol.at(i);
         }

         // Shift old solution storage by one
         for(int n = TimeSchemeSelector::FIELD_MEMORY-1; n > 0; n--)
         {
            rOldSol.at(n) = rOldSol.at(n-1); 
         }
         rOldSol.at(0) = sol;

         // Compute RHS part from nonlinear terms at t_(n-i), i > 0
         for(int n = 0; n < TimeSchemeSelector::nonlinearMemory(step); n++)
         {
            rRHS += TimeSchemeSelector::rhsN(n+1, step)*rOldNL.at(n);
         }

         // Shift old storage by one
         for(int n = TimeSchemeSelector::NONLINEAR_MEMORY-1; n > 0; n--)
         {
            rOldNL.at(n) = rOldNL.at(n-1); 
         }
         rOldNL.at(0) = rTmp;
      }

      template <> inline void computeMVPAY<SparseMatrix,DecoupledZMatrix>(DecoupledZMatrix& y, const SparseMatrix& A, const DecoupledZMatrix& x, const MHDFloat a)
      {
         computeMVPAY<SparseMatrix,Matrix>(y.real(), A, x.real(), a);

         computeMVPAY<SparseMatrix,Matrix>(y.imag(), A, x.imag(), a);
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
         computeMVPAY<SparseMatrix,Matrix>(rRHS.real(), mat, sol.real(), TimeSchemeSelector::rhsN(0, step));

         // Compute RHS part from nonlinear terms at t_(n-i), i > 0
         for(int n = 0; n < TimeSchemeSelector::nonlinearMemory(step); n++)
         {
            rRHS.real() += TimeSchemeSelector::rhsN(n+1, step)*rOldNL.at(n).real();
         }

         // Shift old storage by one
         for(int n = TimeSchemeSelector::NONLINEAR_MEMORY-1; n > 0; n--)
         {
            rOldNL.at(n).real() = rOldNL.at(n-1).real(); 
         }
         rOldNL.at(0).real() = rTmp.real();

         // Imaginary part
         
         // Store RHS at t_n
         rTmp.real() = rRHS.imag();

         // Compute RHS part from solution and nonlinear term at t_n
         computeMVPAY<SparseMatrix,Matrix>(rRHS.imag(), mat, sol.imag(), TimeSchemeSelector::rhsN(0, step));

         // Compute RHS part from nonlinear terms at t_(n-i), i > 0
         for(int n = 0; n < TimeSchemeSelector::nonlinearMemory(step); n++)
         {
            rRHS.imag() += TimeSchemeSelector::rhsN(n+1, step)*rOldNL.at(n).imag();
         }

         // Shift old storage by one
         for(int n = TimeSchemeSelector::NONLINEAR_MEMORY-1; n > 0; n--)
         {
            rOldNL.at(n).imag() = rOldNL.at(n-1).imag(); 
         }
         rOldNL.at(0).imag() = rTmp.real();
      }

      template <> inline void computeRHSNoNLMemory<SparseMatrix,DecoupledZMatrix>(const int step, DecoupledZMatrix& rRHS, const SparseMatrix& mat, const DecoupledZMatrix& sol, const std::vector<SparseMatrix>& oldMat, std::vector<DecoupledZMatrix>& rOldSol)
      {
         // Compute RHS part from solution and nonlinear term at t_n
         computeMVPAY<SparseMatrix,DecoupledZMatrix>(rRHS, mat, sol, TimeSchemeSelector::rhsN(0, step));

         // Real part

         // Compute RHS part from solution at t_(n-i), i > 0
         for(int i = 0; i < TimeSchemeSelector::fieldMemory(step); i++)
         {
            rRHS.real() += oldMat.at(i)*rOldSol.at(i).real();
         }

         // Shift old solution storage by one
         for(int n = TimeSchemeSelector::FIELD_MEMORY-1; n > 0; n--)
         {
            rOldSol.at(n).real() = rOldSol.at(n-1).real(); 
         }
         rOldSol.at(0).real() = sol.real();

         // Imaginary part

         // Compute RHS part from solution at t_(n-i), i > 0
         for(int i = 0; i < TimeSchemeSelector::fieldMemory(step); i++)
         {
            rRHS.imag() += oldMat.at(i)*rOldSol.at(i).imag();
         }

         // Shift old solution storage by one
         for(int n = TimeSchemeSelector::FIELD_MEMORY-1; n > 0; n--)
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
         computeMVPAY<SparseMatrix,Matrix>(rRHS.real(), mat, sol.real(), TimeSchemeSelector::rhsN(0, step));

         // Compute RHS part from solution at t_(n-i), i > 0
         for(int i = 0; i < TimeSchemeSelector::fieldMemory(step); i++)
         {
            rRHS.real() += oldMat.at(i)*rOldSol.at(i).real();
         }

         // Shift old solution storage by one
         for(int n = TimeSchemeSelector::FIELD_MEMORY-1; n > 0; n--)
         {
            rOldSol.at(n).real() = rOldSol.at(n-1).real(); 
         }
         rOldSol.at(0).real() = sol.real();

         // Compute RHS part from nonlinear terms at t_(n-i), i > 0
         for(int n = 0; n < TimeSchemeSelector::nonlinearMemory(step); n++)
         {
            rRHS.real() += TimeSchemeSelector::rhsN(n+1, step)*rOldNL.at(n).real();
         }

         // Shift old storage by one
         for(int n = TimeSchemeSelector::NONLINEAR_MEMORY-1; n > 0; n--)
         {
            rOldNL.at(n).real() = rOldNL.at(n-1).real(); 
         }
         rOldNL.at(0).real() = rTmp.real();

         // Imaginary part
         
         // Store RHS at t_n
         rTmp.real() = rRHS.imag();

         // Compute RHS part from solution and nonlinear term at t_n
         computeMVPAY<SparseMatrix,Matrix>(rRHS.imag(), mat, sol.imag(), TimeSchemeSelector::rhsN(0, step));

         // Compute RHS part from solution at t_(n-i), i > 0
         for(int i = 0; i < TimeSchemeSelector::fieldMemory(step); i++)
         {
            rRHS.imag() += oldMat.at(i)*rOldSol.at(i).imag();
         }

         // Shift old solution storage by one
         for(int n = TimeSchemeSelector::FIELD_MEMORY-1; n > 0; n--)
         {
            rOldSol.at(n).imag() = rOldSol.at(n-1).imag(); 
         }
         rOldSol.at(0).imag() = sol.imag();

         // Compute RHS part from nonlinear terms at t_(n-i), i > 0
         for(int n = 0; n < TimeSchemeSelector::nonlinearMemory(step); n++)
         {
            rRHS.imag() += TimeSchemeSelector::rhsN(n+1, step)*rOldNL.at(n).imag();
         }

         // Shift old storage by one
         for(int n = TimeSchemeSelector::NONLINEAR_MEMORY-1; n > 0; n--)
         {
            rOldNL.at(n).imag() = rOldNL.at(n-1).imag(); 
         }
         rOldNL.at(0).imag() = rTmp.real();
      }
   }
}
}

#endif // SPARSEOLDIMEXTIMESTEPPER_HPP
