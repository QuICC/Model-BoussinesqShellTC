/** \file TimestepCoordinator.hpp
 *  \brief Implementation of a general timestep coordinator structure
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
#include "SparseSolvers/SparseLinearCoordinatorBase.hpp"
#include "Timesteppers/SparseDTimestepper.hpp"
#include "Timesteppers/SparseZTimestepper.hpp"
#include "Timesteppers/ImExRK3.hpp"
#include "Equations/IScalarEquation.hpp"
#include "Equations/IVectorEquation.hpp"

namespace GeoMHDiSCC {

namespace Timestep {

   /**
    * @brief Implementation of general timestepper structure
    */
   class TimestepCoordinator: public Solver::SparseLinearCoordinatorBase
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
          * @brief Create a real linear solver
          */
         virtual void addSolverD(const int start);

         /**
          * @brief Create a complex linear solver
          */
         virtual void addSolverZ(const int start);

         /**
          * @brief Build the real solver matrix
          *
          * @param spSolver   Shared sparse real solver
          * @param matIdx     Index of the solver matrix
          * @param spEq       Shared pointer to equation
          * @param comp       Field component
          * @param idx        Matrix index
          */
         virtual void buildSolverMatrix(Solver::SharedSparseDLinearSolver spSolver, const int matIdx, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const int idx);

         /**
          * @brief Build the complex solver matrix
          *
          * @param spSolver   Shared sparse real solver
          * @param matIdx     Index of the solver matrix
          * @param spEq       Shared pointer to equation
          * @param comp       Field component
          * @param idx        Matrix index
          */
         virtual void buildSolverMatrix(Solver::SharedSparseZLinearSolver spSolver, const int matIdx, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const int idx);

      private:
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
}
}

#endif // TIMESTEPCOORDINATOR_HPP
