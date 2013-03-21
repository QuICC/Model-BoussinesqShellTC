/** \file EquationZTimestepper.hpp
 *  \brief Implementation of a complex valued (coupled) equation timestepper
 *
 *  \mhdBug Needs test
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
#include "Timesteppers/SparseSolverMacro.h"
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
          * @param nField  Number of fields
          * @param start   Starting index (for example without m=0)
          */
         EquationZTimestepper(const int nField, const int start);

         /**
          * @brief Destructor
          */
         virtual ~EquationZTimestepper();

         /**
          * @brief Reserve storage for the matrices
          */
         void reserveMatrices(const int n);

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
          * @brief Add LHS Matrix
          *
          * @param lhs LHS matrix to add
          */
         void addLHSMatrix(const DecoupledZSparse& lhs);

         /**
          * @brief Add RHS Matrix
          *
          * @param rhs RHS matrix to add
          */
         void addRHSMatrix(const DecoupledZSparse& rhs);

         /**
          * @brief Set existing LHS Matrix with new values
          *
          * @param idx  Index of the matrix
          * @param lhs LHS matrix to add
          */
         void setLHSMatrix(const int idx, const DecoupledZSparse& lhs);

         /**
          * @brief Se existing RHS Matrix with new values
          *
          * @param idx  Index of the matrix
          * @param rhs RHS matrix to add
          */
         void setRHSMatrix(const int idx, const DecoupledZSparse& rhs);

         /**
          * @brief Complete LHS Matrix with additional values
          *
          * @param idx  Index of the matrix
          * @param lhs  LHS matrix to add
          */
         void completeLHSMatrix(const int idx, const DecoupledZSparse& lhs);

         /**
          * @brief Complete RHS Matrix with additional values
          *
          * @param idx  Index of the matrix
          * @param rhs RHS matrix to add
          */
         void completeRHSMatrix(const int idx, const DecoupledZSparse& rhs);

         /**
          * @brief Add RHS and solution data storage
          * 
          * @param cols Number of columns required
          */
         void addStorage(const int cols);

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
