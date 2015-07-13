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
#include "Enums/ModelOperator.hpp"
#include "Enums/ModelOperatorBoundary.hpp"
#include "SparseSolvers/SparseLinearCoordinatorBase.hpp"
#include "Timesteppers/SparseTimestepper.hpp"
#include "Equations/IScalarEquation.hpp"
#include "Equations/IVectorEquation.hpp"
#include "IoAscii/CflWriter.hpp"

namespace GeoMHDiSCC {

namespace Timestep {

   /**
    * @brief Implementation of general timestepper structure
    */
   class TimestepCoordinator: public Solver::SparseLinearCoordinatorBase<SparseTimestepper>
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
          * @param time    Initial time value
          * @param dt      Initial timestep value
          * @param scalEq  Shared scalar equations
          * @param vectEq  Shared vector equations
          */
         void init(const MHDFloat time, const MHDFloat dt, const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq);

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
          * @brief Build the real operator, complex field solver matrix
          *
          * @param spSolver   Shared sparse real solver
          * @param matIdx     Index of the solver matrix
          * @param spEq       Shared pointer to equation
          * @param comp       Field component
          * @param idx        Matrix index
          */
         virtual void buildSolverMatrix(TimestepCoordinator::SharedRealSolverType spSolver, const int matIdx, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const int idx);

         /**
          * @brief Build the complex operator, complex field solver matrix
          *
          * @param spSolver   Shared sparse real solver
          * @param matIdx     Index of the solver matrix
          * @param spEq       Shared pointer to equation
          * @param comp       Field component
          * @param idx        Matrix index
          */
         virtual void buildSolverMatrix(TimestepCoordinator::SharedComplexSolverType spSolver, const int matIdx, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const int idx);

      private:
         using Solver::SparseLinearCoordinatorBase<SparseTimestepper>::init;

         /**
          * @brief Compute the timestep operator coefficients
          */
         void computeTimeCoeffs(MHDFloat& lhsL, MHDFloat& lhsT, MHDFloat& rhsL, MHDFloat& rhsT, std::vector<MHDFloat>& oldRhsL, std::vector<MHDFloat>& oldRhsT);

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
          * @brief No update window for timestep decrease
          */
         const MHDFloat mcDownWindow;

         /**
          * @brief Minimal timestep allowed before simulation abort
          */
         const MHDFloat mcMinDt;

         /**
          * @brief Maximum timestep allowed 
          */
         const MHDFloat mcMaxDt;

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

         /**
          * @brief Shared CFL writer
          */
         IoAscii::SharedCflWriter   mspIo;
   };

   /**
    * @brief Small wrapper for a generic implementation of the solver matrix construction
    */
   template <typename TStepper> void buildSolverMatrixWrapper(typename SharedPtrMacro<TStepper > spSolver, const int matIdx, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const int idx, const MHDFloat lhsLCoeff, const MHDFloat lhsTCoeff, const MHDFloat rhsLCoeff, const MHDFloat rhsTCoeff, const std::vector<MHDFloat>& oldRhsLCoeff, const std::vector<MHDFloat>& oldRhsTCoeff);

   template <typename TStepper> void buildSolverMatrixWrapper(typename SharedPtrMacro<TStepper > spSolver, const int matIdx, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const int idx, const MHDFloat lhsLCoeff, const MHDFloat lhsTCoeff, const MHDFloat rhsLCoeff, const MHDFloat rhsTCoeff, const std::vector<MHDFloat>& oldRhsLCoeff, const std::vector<MHDFloat>& oldRhsTCoeff)
   {
      // Resize LHS matrix if necessary
      spSolver->rLHSMatrix(matIdx).resize(spEq->couplingInfo(comp).systemN(idx), spEq->couplingInfo(comp).systemN(idx));

      // Resize RHS matrix at t_n if necessary
      spSolver->rRHSMatrix(matIdx).resize(spEq->couplingInfo(comp).systemN(idx), spEq->couplingInfo(comp).systemN(idx));

      // Resize RHS matrix at t_(n-i), i > 0 if necessary
      for(size_t i = 0; i < oldRhsLCoeff.size(); i++)
      {
         spSolver->rOldRHSMatrix(i, matIdx).resize(spEq->couplingInfo(comp).systemN(idx), spEq->couplingInfo(comp).systemN(idx));
      }

      DecoupledZSparse  linOp;
      DecoupledZSparse  timeOp;

      // Compute model's linear operator
      spEq->buildModelMatrix(linOp, ModelOperator::IMPLICIT_LINEAR, comp, idx, ModelOperatorBoundary::SOLVER_HAS_BC);
      // Compute model's time operator
      spEq->buildModelMatrix(timeOp, ModelOperator::TIME, comp, idx, ModelOperatorBoundary::SOLVER_NO_TAU);

      // Set LHS matrix
      Solver::internal::addOperators(spSolver->rLHSMatrix(matIdx), lhsLCoeff, linOp); 
      Solver::internal::addOperators(spSolver->rLHSMatrix(matIdx), -lhsTCoeff, timeOp); 

      // Set RHS matrix at t_n
      spEq->buildModelMatrix(linOp, ModelOperator::IMPLICIT_LINEAR, comp, idx, ModelOperatorBoundary::SOLVER_NO_TAU);
      spEq->buildModelMatrix(timeOp, ModelOperator::TIME, comp, idx, ModelOperatorBoundary::SOLVER_NO_TAU);
      Solver::internal::addOperators(spSolver->rRHSMatrix(matIdx), -rhsLCoeff, linOp); 
      Solver::internal::addOperators(spSolver->rRHSMatrix(matIdx), -rhsTCoeff, timeOp); 

      // Set RHS matrix at t_(n-i), i > 0
      for(size_t i = 0; i < oldRhsLCoeff.size(); i++)
      {
         Solver::internal::addOperators(spSolver->rOldRHSMatrix(i, matIdx), -oldRhsLCoeff.at(i), linOp); 
         Solver::internal::addOperators(spSolver->rOldRHSMatrix(i, matIdx), -oldRhsTCoeff.at(i), timeOp); 
      }

      // Set time matrix for timestep updates
      if(matIdx == idx)
      {
         // Resize time matrix if necessary
         spSolver->rTMatrix(idx).resize(spEq->couplingInfo(comp).systemN(idx), spEq->couplingInfo(comp).systemN(idx));

         // Set time matrix
         Solver::internal::addOperators(spSolver->rTMatrix(idx), 1.0, timeOp); 
      }
      
      // Solver is initialized
      spSolver->setInitialized();
   }

   // 
   // Dummy solver specialization
   //

   template <> inline void buildSolverMatrixWrapper<Solver::SparseDummySolverComplexType>(SharedPtrMacro<Solver::SparseDummySolverComplexType> spSolver, const int matIdx, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const int idx, const MHDFloat lhsLCoeff, const MHDFloat lhsTCoeff, const MHDFloat rhsLCoeff, const MHDFloat rhsTCoeff, const std::vector<MHDFloat>& oldRhsLCoeff, const std::vector<MHDFloat>& oldRhsTCoeff) {};

}
}

#endif // TIMESTEPCOORDINATOR_HPP
