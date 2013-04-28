/** \file EquationZTimestepper.hpp
 *  \brief Implementation of a complex valued (coupled) equation timestepper
 */

#ifndef EQUATIONZTIMESTEPPER_HPP
#define EQUATIONZTIMESTEPPER_HPP

// Configuration includes
//
#include "SmartPointers/SharedPtrMacro.h"

// System includes
//

// External includes
//

// Project includes
//
#include "SparseSolvers/SparseSolverMacro.h"
#include "Timesteppers/EquationTimestepperBase.hpp"

namespace GeoMHDiSCC {

namespace Timestep {

   /**
    *  \brief Implementation of a complex valued (coupled) equation timestepper
    */
   class EquationZTimestepper: public EquationTimestepperBase
   {
      public:
         /**
          * @brief Constructor
          *
          * @param start   Starting index (for example without m=0)
          */
         EquationZTimestepper(const int start);

         /**
          * @brief Destructor
          */
         virtual ~EquationZTimestepper();

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
          * @brief Update the LHS matrix with new timedependence
          *
          * @param lhsCoeff   New coefficient for LHS time dependent part
          * @param rhsCoeff   New coefficient for RHS time dependent part
          * @param step       Timestep scheme substep
          */
         void updateTimeMatrix(const MHDFloat lhsCoeff, const MHDFloat rhsCoeff, const int step);

         /**
          * @brief Get the number of linear systems in solver
          */
         int nSystem() const;

         /**
          * @brief Set LHS matrix
          *
          * @param idx Index of the matrix
          */
         SparseMatrixZ& rLHSMatrix(const int idx);

         /**
          * @brief Set RHS matrix
          *
          * @param idx Index of the matrix
          */
         SparseMatrixZ& rRHSMatrix(const int idx);

         /**
          * @brief Set time dependent part of LHS matrix
          *
          * @param idx Index of the matrix
          */
         SparseMatrixZ& rTMatrix(const int idx);

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
         MatrixZ& rRHSData(const int idx);

         /**
          * @brief Get solution data
          *
          * @param idx   Index of the data
          */
         const MatrixZ& solution(const int idx) const;

         /**
          * @brief Set solution data
          *
          * @param idx   Index of the data
          */
         MatrixZ& rSolution(const int idx);
         
      protected:
         /**
          * @brief Complex LHS operators of the timestepped equations
          */
         std::vector<SparseMatrixZ>   mLHSMatrix;

         /**
          * @brief Complex RHS operators of the timestepped equations
          */
         std::vector<SparseMatrixZ>   mRHSMatrix;

         /**
          * @brief Time dependent part of LHS matrix
          */
         std::vector<SparseMatrixZ>   mTMatrix;

         /**
          * @brief Storage for linear solve's RHS
          */
         std::vector<MatrixZ>  mRHSData;

         /**
          * @brief Storage for old nonlinear RHS
          */
         std::vector<MatrixZ>  mRHSOld;

         /**
          * @brief Storage for solution of linear solve
          */
         std::vector<MatrixZ>  mSolution;

         /**
          * @brief Create sparse solvers
          */
         std::vector<SharedPtrMacro<SparseSolverMacro<SparseMatrixZ> > >  mSolver;

      private:
   };

   inline int EquationZTimestepper::nSystem() const
   {
      return this->mRHSData.size();
   }

   inline MatrixZ& EquationZTimestepper::rRHSData(const int idx)
   {
      return this->mRHSData.at(idx);
   }

   inline const MatrixZ& EquationZTimestepper::solution(const int idx) const
   {
      return this->mSolution.at(idx);
   }

   inline MatrixZ& EquationZTimestepper::rSolution(const int idx)
   {
      return this->mSolution.at(idx);
   }

   /// Typedef for a shared pointer of a EquationZTimestepper
   typedef SharedPtrMacro<EquationZTimestepper>  SharedEquationZTimestepper;
}
}

#endif // EQUATIONZTIMESTEPPER_HPP
