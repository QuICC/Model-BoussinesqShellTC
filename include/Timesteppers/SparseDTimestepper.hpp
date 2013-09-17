/** 
 * @file SparseDTimestepper.hpp
 * @brief Implementation of a real valued (coupled) equation timestepper
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef SPARSEDTIMESTEPPER_HPP
#define SPARSEDTIMESTEPPER_HPP

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
#include "SparseSolvers/SparseDLinearSolver.hpp"

namespace GeoMHDiSCC {

namespace Timestep {

   /**
    * @brief Implementation of a real valued (coupled) equation timestepper
    */
   class SparseDTimestepper: public Solver::SparseDLinearSolver
   {
      public:
         /**
          * @brief Constructor
          *
          * @param start   Starting index (for example without m=0)
          */
         SparseDTimestepper(const int start);

         /**
          * @brief Destructor
          */
         virtual ~SparseDTimestepper();

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
          * @brief Update the LHS matrix triplets
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
         virtual void addStorage(const int rows, const int cols);
         
      protected:
         /**
          * @brief Real RHS operators of the timestepped equations
          */
         std::vector<SparseMatrix>   mRHSMatrix;

         /**
          * @brief Time dependent part of LHS matrix
          */
         std::vector<SparseMatrix>   mTMatrix;

         /**
          * @brief Storage for old nonlinear RHS
          */
         std::vector<DecoupledZMatrix>  mRHSOld;

      private:
   };

   /// Typedef for a shared pointer of a SparseDTimestepper
   typedef SharedPtrMacro<SparseDTimestepper>  SharedSparseDTimestepper;
}
}

#endif // SPARSEDTIMESTEPPER_HPP
