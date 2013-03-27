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
          * @brief Set LHS triplets
          *
          * @param idx Index of the triplets
          */
         std::vector<TripletZ>& rLHSTriplets(const int idx);

         /**
          * @brief Set RHS triplets
          *
          * @param idx Index of the triplets
          */
         std::vector<TripletZ>& rRHSTriplets(const int idx);

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
          * @brief Storage for the matrix sizes
          */
         std::vector<int>  mSize;

         /**
          * @brief Complex RHS operators of the timestepped equations
          */
         std::vector<SparseMatrixZ>   mRHSMatrix;

         /**
          * @brief Complex LHS operators triplets of the timestepped equations
          */
         std::vector<std::vector<TripletZ> >   mLHSTriplets;

         /**
          * @brief Complex RHS operators triplets of the timestepped equations
          */
         std::vector<std::vector<TripletZ> >   mRHSTriplets;

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
