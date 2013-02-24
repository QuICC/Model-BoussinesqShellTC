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
#include <Eigen/SuperLUSupport>

// Project includes
//
#include "Timesteppers/EquationTimestepperBase.hpp"

namespace GeoMHDiSCC {

   /**
    *  \brief Implementation of a complex valued (coupled) equation timestepper
    */
   class EquationZTimestepper: public EquationTimestepperBase
   {
      public:
         /**
          * @brief Constructor
          */
         EquationZTimestepper(const int nField);

         /**
          * @brief Destructor
          */
         virtual ~EquationZTimestepper();

         /**
          * @brief Initialise solvers
          */
         void initSolver();

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
          * @para idx   Index of the data
          */
         MatrixZ& rRHSData(const int i);

         /**
          * @brief Get solution data
          *
          * @para idx   Index of the data
          */
         const MatrixZ& solution(const int i) const;

         /**
          * @brief Set solution data
          *
          * @para idx   Index of the data
          */
         MatrixZ& rSolution(const int i);
         
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
         std::vector<SharedPtrMacro<Eigen::SuperLU<SparseMatrixZ> > >  mSolver;

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

#endif // EQUATIONZTIMESTEPPER_HPP
