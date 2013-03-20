/** \file EquationDTimestepper.hpp
 *  \brief Implementation of a real valued (coupled) equation timestepper
 *
 *  \mhdBug Needs test
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
          */
         EquationDTimestepper(const int nField);

         /**
          * @brief Destructor
          */
         virtual ~EquationDTimestepper();

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
          * @param lhs LHS matrix to add
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
