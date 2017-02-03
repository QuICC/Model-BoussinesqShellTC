/** 
 * @file SparseLinearCoordinator.hpp
 * @brief Implementation of a general sparse linear solver coordinator
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef SPARSELINEARCOORDINATOR_HPP
#define SPARSELINEARCOORDINATOR_HPP

// Configuration includes
//
#include "SmartPointers/SharedPtrMacro.h"

// System includes
//

// External includes
//

// Project includes
//
#include "Enums/ModelOperator.hpp"
#include "Enums/ModelOperatorBoundary.hpp"
#include "SparseSolvers/SparseLinearCoordinatorBase.hpp"
#include "SparseSolvers/SparseLinearSolver.hpp"
#include "Equations/IScalarEquation.hpp"
#include "Equations/IVectorEquation.hpp"

namespace QuICC {

namespace Solver {

   /**
    * @brief Implementation of general sparse linear solver coordinator
    */
   class SparseLinearCoordinator: public SparseLinearCoordinatorBase<SparseLinearSolver>
   {
      public:
         /**
          * @brief Constructor
          */
         SparseLinearCoordinator();

         /**
          * @brief Destructor
          */
         ~SparseLinearCoordinator();

         /**
          * @brief Solve the equations
          *
          * @param scalEq Shared scalar equations
          * @param vectEq Shared vector equations
          * @param scalVar Shared scalar variables
          * @param vectVar Shared vector variables
          */
         void solve(const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq, const ScalarVariable_map& scalVar, const VectorVariable_map& vectVar);

      protected:
         /**
          * @brief Build the real operator, complex field solver matrix
          *
          * @param spSolver   Shared sparse real solver
          * @param spEq       Shared pointer to equation
          * @param comp       Field component
          * @param idx        Matrix index
          */
         virtual void buildSolverMatrix(SparseLinearCoordinator::SharedRealSolverType spSolver, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const int idx);

         /**
          * @brief Build the complex operator, complex field solver matrix
          *
          * @param spSolver   Shared sparse real solver
          * @param spEq       Shared pointer to equation
          * @param comp       Field component
          * @param idx        Matrix index
          */
         virtual void buildSolverMatrix(SparseLinearCoordinator::SharedComplexSolverType spSolver, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const int idx);

      private:
   };

   /**
    * @brief Generic implementation to build solver matrices
    */
   template <typename TSolver> void buildLinearSolverMatrixWrapper(SharedPtrMacro<TSolver > spSolver, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const int idx);

   template <typename TSolver> void buildLinearSolverMatrixWrapper(SharedPtrMacro<TSolver > spSolver, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const int idx)
   {
      DecoupledZSparse  linOp;

      // Get model's linear operator with tau lines
      spEq->buildModelMatrix(linOp, ModelOperator::IMPLICIT_LINEAR, comp, idx, ModelOperatorBoundary::SOLVER_HAS_BC);

      spSolver->buildOperators(idx, linOp, spEq->couplingInfo(comp).systemN(idx));

      // Solver is initialized
      spSolver->setInitialized();
   }

   //
   // Dummy solver specialization
   //
   
   template <> inline void buildLinearSolverMatrixWrapper<SparseDummySolverComplexType>(SharedPtrMacro<SparseDummySolverComplexType > spSolver, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const int idx) {};
}
}

#endif // SPARSELINEARCOORDINATOR_HPP
