/** \file SparseZTimestepper.hpp
 *  \brief Implementation of a complex valued (coupled) equation timestepper
 */

#ifndef SPARSEZTIMESTEPPER_HPP
#define SPARSEZTIMESTEPPER_HPP

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
#include "SparseSolvers/SparseZLinearSolver.hpp"

namespace GeoMHDiSCC {

namespace Timestep {

   /**
    *  \brief Implementation of a complex valued (coupled) equation timestepper
    */
   class SparseZTimestepper: public Solver::SparseZLinearSolver
   {
      public:
         /**
          * @brief Constructor
          *
          * @param start   Starting index (for example without m=0)
          */
         SparseZTimestepper(const int start);

         /**
          * @brief Destructor
          */
         virtual ~SparseZTimestepper();

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
          * @param lhsCoeff   New coefficient for LHS time dependent part
          * @param rhsCoeff   New coefficient for RHS time dependent part
          * @param step       Timestep scheme substep
          */
         void updateTimeMatrix(const MHDFloat lhsCoeff, const MHDFloat rhsCoeff, const int step);

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
         virtual void addStorage(const int rows, const int cols);
         
      protected:
         /**
          * @brief Complex RHS operators of the timestepped equations
          */
         std::vector<SparseMatrixZ>   mRHSMatrix;

         /**
          * @brief Time dependent part of LHS matrix
          */
         std::vector<SparseMatrixZ>   mTMatrix;

         /**
          * @brief Storage for old nonlinear RHS
          */
         std::vector<MatrixZ>  mRHSOld;

      private:
   };

   /// Typedef for a shared pointer of a SparseZTimestepper
   typedef SharedPtrMacro<SparseZTimestepper>  SharedSparseZTimestepper;
}
}

#endif // SPARSEZTIMESTEPPER_HPP
