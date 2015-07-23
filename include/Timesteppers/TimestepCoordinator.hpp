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
#include "TypeSelectors/TimeSchemeTypeSelector.hpp"
#include "Equations/IScalarEquation.hpp"
#include "Equations/IVectorEquation.hpp"
#include "IoAscii/CflWriter.hpp"

namespace GeoMHDiSCC {

namespace Timestep {

   /**
    * @brief Implementation of general timestepper structure
    */
   class TimestepCoordinator: public Solver::SparseLinearCoordinatorBase<TimeSchemeTypeSelector>
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
          * @brief Tune adaptive timestepper 
          *
          * \mhdBug Not fully implemented
          */
         void tuneAdaptive(const MHDFloat time);

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
          * @param spEq       Shared pointer to equation
          * @param comp       Field component
          * @param idx        Matrix index
          */
         virtual void buildSolverMatrix(TimestepCoordinator::SharedRealSolverType spSolver, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const int idx);

         /**
          * @brief Build the complex operator, complex field solver matrix
          *
          * @param spSolver   Shared sparse real solver
          * @param spEq       Shared pointer to equation
          * @param comp       Field component
          * @param idx        Matrix index
          */
         virtual void buildSolverMatrix(TimestepCoordinator::SharedComplexSolverType spSolver, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const int idx);

      private:
         using Solver::SparseLinearCoordinatorBase<TimeSchemeTypeSelector>::init;

         /**
          * @brief Update time dependence
          */
         void updateMatrices();

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
          * @brief Constant timestep steps
          */
         MHDFloat mCnstSteps;

         /**
          * @brief Constant timestep steps
          */
         MHDFloat mStepTime;

         /**
          * @brief Shared CFL writer
          */
         IoAscii::SharedCflWriter   mspIo;
   };

   /**
    * @brief Small wrapper for a generic implementation of the solver matrix construction
    */
   template <typename TStepper> void buildTimestepMatrixWrapper(typename SharedPtrMacro<TStepper > spSolver, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const MHDFloat dt, const int idx);

   template <typename TStepper> void buildTimestepMatrixWrapper(typename SharedPtrMacro<TStepper > spSolver, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const MHDFloat dt, const int idx)
   {
      DecoupledZSparse  linOp;
      DecoupledZSparse  timeOp;
      DecoupledZSparse  bcOp;

      // Compute model's linear operator (without Tau lines)
      spEq->buildModelMatrix(linOp, ModelOperator::IMPLICIT_LINEAR, comp, idx, ModelOperatorBoundary::SOLVER_NO_TAU);
      // Compute model's time operator (without Tau lines)
      spEq->buildModelMatrix(timeOp, ModelOperator::TIME, comp, idx, ModelOperatorBoundary::SOLVER_NO_TAU);
      // Compute model's tau line boundary operator
      spEq->buildModelMatrix(bcOp, ModelOperator::BOUNDARY, comp, idx, ModelOperatorBoundary::SOLVER_HAS_BC);

      // Let the timestepper build the right operators
      spSolver->buildOperators(idx, linOp, timeOp, bcOp, dt, spEq->couplingInfo(comp).systemN(idx));
      
      // Solver is initialized
      spSolver->setInitialized();
   }

   // 
   // Dummy solver specialization
   //

   template <> inline void buildTimestepMatrixWrapper<Solver::SparseDummySolverComplexType>(SharedPtrMacro<Solver::SparseDummySolverComplexType> spSolver, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const MHDFloat dt, const int idx) {};

}
}

#endif // TIMESTEPCOORDINATOR_HPP
