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

namespace GeoMHDiSCC {

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
          * @param matIdx     Index of the solver matrix
          * @param spEq       Shared pointer to equation
          * @param comp       Field component
          * @param idx        Matrix index
          */
         virtual void buildSolverMatrix(SparseLinearCoordinator::SharedRealSolverType spSolver, const int matIdx, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const int idx);

         /**
          * @brief Build the complex operator, complex field solver matrix
          *
          * @param spSolver   Shared sparse real solver
          * @param matIdx     Index of the solver matrix
          * @param spEq       Shared pointer to equation
          * @param comp       Field component
          * @param idx        Matrix index
          */
         virtual void buildSolverMatrix(SparseLinearCoordinator::SharedComplexSolverType spSolver, const int matIdx, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const int idx);

      private:
   };

   /**
    * @brief Generic implementation to build solver matrices
    */
   template <typename TSolver> void buildLinearSolverMatrixWrapper(SharedPtrMacro<TSolver > spSolver, const int matIdx, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const int idx);

   template <typename TSolver> void buildLinearSolverMatrixWrapper(SharedPtrMacro<TSolver > spSolver, const int matIdx, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const int idx)
   {
      // Resize operator matrix
      spSolver->rLHSMatrix(matIdx).resize(spEq->couplingInfo(comp).systemN(idx), spEq->couplingInfo(comp).systemN(idx));

      // Get model operator
      DecoupledZSparse  linOp;
      spEq->buildModelMatrix(linOp, ModelOperator::IMPLICIT_LINEAR, comp, idx, ModelOperatorBoundary::HAS_BC);
      Solver::internal::addOperators(spSolver->rLHSMatrix(matIdx), 1.0, linOp); 

      // Solver is initialized
      spSolver->setInitialized();
   }

   //
   // Dummy solver specialization
   //
   
   template <> inline void buildLinearSolverMatrixWrapper<SparseDummySolverComplexType>(SharedPtrMacro<SparseDummySolverComplexType > spSolver, const int matIdx, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const int idx) {};
}
}

#endif // SPARSELINEARCOORDINATOR_HPP
