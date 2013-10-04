/** 
 * @file TimestepCoordinator.hpp
 * @brief Implementation of a general timestep coordinator structure
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef TIMESTEPCOORDINATOR_HPP
#define TIMESTEPCOORDINATOR_HPP

// Configuration includes
//
#include "SmartPointers/SharedPtrMacro.h"

// System includes
//

// External includes
//

// Project includes
//
#include "Base/MathConstants.hpp"
#include "SparseSolvers/SparseLinearCoordinatorBase.hpp"
#include "Timesteppers/SparseTimestepper.hpp"
#include "Equations/IScalarEquation.hpp"
#include "Equations/IVectorEquation.hpp"

namespace GeoMHDiSCC {

namespace Timestep {

   /**
    * @brief Implementation of general timestepper structure
    */
   class TimestepCoordinator: public Solver::SparseLinearCoordinatorBase<SparseZTimestepper,SparseRZTimestepper>
   {
      public:
         /// Typedef for a shared scalar equation iterator
         typedef std::vector<Equations::SharedIScalarEquation>::iterator   ScalarEquation_iterator;

         /// Typedef for a shared vector equation iterator
         typedef std::vector<Equations::SharedIVectorEquation>::iterator   VectorEquation_iterator;

         /// Typedef for a shared scalar equation range
         typedef std::pair<ScalarEquation_iterator, ScalarEquation_iterator>  ScalarEquation_range;

         /// Typedef for a shared vector equation range
         typedef std::pair<VectorEquation_iterator, VectorEquation_iterator>  VectorEquation_range;

         /**
          * @brief Constructor
          */
         TimestepCoordinator();

         /**
          * @brief Destructor
          */
         ~TimestepCoordinator();

         /**
          * @brief Initialise timestepper
          *
          * @param dt      Initial timestep value
          * @param scalEq  Shared scalar equations
          * @param vectEq  Shared vector equations
          */
         void init(const MHDFloat dt, const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq);

         /**
          * @brief Adapt the timestep used
          *
          * @param cfl     CFL condition
          * @param scalEq  Shared scalar equations
          * @param vectEq  Shared vector equations
          */
         void adaptTimestep(const MHDFloat cfl, const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq);

         /**
          * @brief Update control status
          */
         void update();

         /**
          * @brief Compute (partial) forward step
          *
          * @param scalEq Shared scalar equations
          * @param vectEq Shared vector equations
          * @param scalVar Shared scalar variables
          * @param vectVar Shared vector variables
          */
         void stepForward(const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq, const ScalarVariable_map& scalVar, const VectorVariable_map& vectVar);

         /**
          * @brief Get current simulation time
          */
         MHDFloat time() const;

         /**
          * @brief Get current simulation timestep
          */
         MHDFloat timestep() const;
         
      protected:
         /**
          * @brief Build the real solver matrix
          *
          * @param spSolver   Shared sparse real solver
          * @param matIdx     Index of the solver matrix
          * @param spEq       Shared pointer to equation
          * @param comp       Field component
          * @param idx        Matrix index
          */
         virtual void buildSolverMatrix(TimestepCoordinator::SharedRSolverType spSolver, const int matIdx, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const int idx);

         /**
          * @brief Build the complex solver matrix
          *
          * @param spSolver   Shared sparse real solver
          * @param matIdx     Index of the solver matrix
          * @param spEq       Shared pointer to equation
          * @param comp       Field component
          * @param idx        Matrix index
          */
         virtual void buildSolverMatrix(TimestepCoordinator::SharedZSolverType spSolver, const int matIdx, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const int idx);

      private:
         template <typename TStepper> void buildSolverMatrixWrapper(typename SharedPtrMacro<TStepper > spSolver, const int matIdx, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const int idx);

         /**
          * @brief Update time dependence
          */
         void updateMatrices();

         /**
          * @brief Compute the RHS of all linear systems
          */
         void computeRHS();

         /**
          * @brief Maximum timestep jump per step (See Soederlind)
          */
         const MHDFloat mcMaxJump;

         /**
          * @brief No update window for timestep increase
          */
         const MHDFloat mcUpWindow;

         /**
          * @brief Minimal timestep allowed before simulation abort
          */
         const MHDFloat mcMinDt;

         /**
          * @brief Previous timestep length
          */
         MHDFloat mOldDt;

         /**
          * @brief Timestep length
          */
         MHDFloat mDt;

         /**
          * @brief Current time
          */
         MHDFloat mTime;
   };

   template <typename TStepper> void TimestepCoordinator::buildSolverMatrixWrapper(typename SharedPtrMacro<TStepper > spSolver, const int matIdx, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const int idx)
   {
      // Operator coefficients
      MHDFloat lhsTCoeff;
      MHDFloat rhsTCoeff;
      MHDFloat lhsLCoeff;
      MHDFloat rhsLCoeff;

      // Set time coefficients for LHS Matrix
      lhsTCoeff = TimeSchemeType::lhsT(this->mStep)*1.0/this->mDt;

      // Set linear coefficients for LHS Matrix
      lhsLCoeff = TimeSchemeType::lhsL(this->mStep);

      // Set time coefficients for RHS Matrix
      rhsTCoeff = TimeSchemeType::rhsT(this->mStep)*1.0/this->mDt;

      // Set linear coefficients for RHS Matrix
      rhsLCoeff = -TimeSchemeType::rhsL(this->mStep);

      // Resize LHS matrix if necessary
      if(spSolver->rLHSMatrix(matIdx).size() == 0)
      {
         spSolver->rLHSMatrix(matIdx).resize(spEq->couplingInfo(comp).systemN(idx), spEq->couplingInfo(comp).systemN(idx));
      }

      // Resize RHS matrix if necessary
      if(spSolver->rRHSMatrix(matIdx).size() == 0)
      {
         spSolver->rRHSMatrix(matIdx).resize(spEq->couplingInfo(comp).systemN(idx), spEq->couplingInfo(comp).systemN(idx));
      }

      // Add boundary row for LHS operator
      Solver::internal::addRow(spSolver->rLHSMatrix(matIdx), 1.0, spEq->operatorRow(Equations::IEquation::BOUNDARYROW, comp, idx));

      DecoupledZSparse linRow = spEq->operatorRow(Equations::IEquation::LINEARROW, comp, idx);
      DecoupledZSparse tRow = spEq->operatorRow(Equations::IEquation::TIMEROW, comp, idx);

      // Set LHS matrix
      Solver::internal::addRow(spSolver->rLHSMatrix(matIdx), lhsLCoeff, linRow);
      Solver::internal::addRow(spSolver->rLHSMatrix(matIdx), -lhsTCoeff, tRow);

      // Set RHS matrix
      Solver::internal::addRow(spSolver->rRHSMatrix(matIdx), rhsLCoeff, linRow);
      Solver::internal::addRow(spSolver->rRHSMatrix(matIdx), -rhsTCoeff, tRow);

      // Set time matrix for timestep updates
      if(matIdx == idx)
      {
         // Resize time matrix if necessary
         if(spSolver->rTMatrix(idx).size() == 0)
         {
            spSolver->rTMatrix(idx).resize(spEq->couplingInfo(comp).systemN(idx), spEq->couplingInfo(comp).systemN(idx));
         }

         // Set time matrix
         Solver::internal::addRow(spSolver->rTMatrix(idx), 1.0, tRow);
      }
   }
}
}

#endif // TIMESTEPCOORDINATOR_HPP
