/** \file EquationDTimestepper.hpp
 *  \brief Implementation of a real valued (coupled) equation timestepper
 */

#ifndef EQUATIONDTIMESTEPPER_HPP
#define EQUATIONDTIMESTEPPER_HPP

// Configuration includes
//
#include "SmartPointers/SharedPtrMacro.h"

// System includes
//

// External includes
//

// Project includes
//
#include "Timesteppers/SparseSolverMacro.h"
#include "Timesteppers/EquationTimestepperBase.hpp"

namespace GeoMHDiSCC {

namespace Timestep {

   /**
    *  \brief Implementation of a real valued (coupled) equation timestepper
    */
   class EquationDTimestepper: public EquationTimestepperBase
   {
      public:
         /**
          * @brief Constructor
          *
          * @param start   Starting index (for example without m=0)
          */
         EquationDTimestepper(const int start);

         /**
          * @brief Destructor
          */
         virtual ~EquationDTimestepper();

         /**
          * @brief Initialise the solver matrices storage
          *
          * @param n Size of matrices
          */
         void initMatrices(const int n);

         /**
          * @brief Initialise solver
          */
         void initSolver();

         /**
          * @brief Update solver
          */
         void updateSolver();

         /**
          * @brief Compute the RHS of the linear systems
          *
          * @param step    Current substep
          */
         void computeRHS(const int step);

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
         SparseMatrix& rLHSMatrix(const int idx);

         /**
          * @brief Set RHS matrix
          *
          * @param idx Index of the matrix
          */
         SparseMatrix& rRHSMatrix(const int idx);

         /**
          * @brief Set time dependent part of LHS matrix
          *
          * @param idx Index of the matrix
          */
         SparseMatrix& rTMatrix(const int idx);

         /**
          * @brief Add RHS and solution data storage
          * 
          * @param rows Number of rows of matrix
          * @param cols Number of columns required
          */
         void addStorage(const int rows, const int cols);

         /**
          * @brief Set RHS data
          *
          * @param idx   Index of the data
          */
         DecoupledZMatrix& rRHSData(const int idx);

         /**
          * @brief Get solution data
          *
          * @param idx   Index of the data
          */
         const DecoupledZMatrix& solution(const int idx) const;

         /**
          * @brief Set solution data
          *
          * @param idx   Index of the data
          */
         DecoupledZMatrix& rSolution(const int idx);
         
      protected:
         /**
          * @brief Real LHS operators of the timestepped equations
          */
         std::vector<SparseMatrix>   mLHSMatrix;

         /**
          * @brief Real RHS operators of the timestepped equations
          */
         std::vector<SparseMatrix>   mRHSMatrix;

         /**
          * @brief Time dependent part of LHS matrix
          */
         std::vector<SparseMatrix>   mTMatrix;

         /**
          * @brief Storage for linear solve's RHS
          */
         std::vector<DecoupledZMatrix>  mRHSData;

         /**
          * @brief Storage for old nonlinear RHS
          */
         std::vector<DecoupledZMatrix>  mRHSOld;

         /**
          * @brief Storage for solution of linear solve
          */
         std::vector<DecoupledZMatrix>  mSolution;

         /**
          * @brief Create sparse solvers
          */
         std::vector<SharedPtrMacro<SparseSolverMacro<SparseMatrix> > >  mSolver;

      private:
   };

   inline int EquationDTimestepper::nSystem() const
   {
      return this->mRHSData.size();
   }

   inline DecoupledZMatrix& EquationDTimestepper::rRHSData(const int idx)
   {
      return this->mRHSData.at(idx);
   }

   inline const DecoupledZMatrix& EquationDTimestepper::solution(const int idx) const
   {
      return this->mSolution.at(idx);
   }

   inline DecoupledZMatrix& EquationDTimestepper::rSolution(const int idx)
   {
      return this->mSolution.at(idx);
   }

   /// Typedef for a shared pointer of a EquationDTimestepper
   typedef SharedPtrMacro<EquationDTimestepper>  SharedEquationDTimestepper;
}
}

#endif // EQUATIONDTIMESTEPPER_HPP
