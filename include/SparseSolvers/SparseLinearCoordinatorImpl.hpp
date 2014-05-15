/** 
 * @file SparseLinearCoordinatorImpl.hpp
 * @brief Implementation of a general sparse linear solver coordinator
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef SPARSELINEARCOORDINATORIMPL_HPP
#define SPARSELINEARCOORDINATORIMPL_HPP

// Configuration includes
//
#include "SmartPointers/SharedPtrMacro.h"

// System includes
//

// External includes
//

// Project includes
//
#include "TypeSelectors/ScalarSelector.hpp"
#include "Enums/ModelOperator.hpp"
#include "SparseSolvers/SparseLinearCoordinatorBase.hpp"
#include "SparseSolvers/SparseLinearSolver.hpp"
#include "Equations/IScalarEquation.hpp"
#include "Equations/IVectorEquation.hpp"

namespace GeoMHDiSCC {

namespace Solver {

   /**
    * @brief Implementation of general sparse linear solver coordinator
    */
   class SparseLinearCoordinatorImpl: public SparseLinearCoordinatorBase<SparseLinearSolver,Datatypes::ScalarSelector<Dimensions::Transform::TRA1D>::FIELD_IS_COMPLEX>
   {
      public:
         /**
          * @brief Constructor
          */
         SparseLinearCoordinatorImpl();

         /**
          * @brief Destructor
          */
         ~SparseLinearCoordinatorImpl();

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

      private:
         /**
          * @brief Wrapper to allow for a unique implementation of all the buildSolverMatrix methods
          */
         template <typename TStepper> void buildSolverMatrixWrapper(SharedPtrMacro<TStepper > spSolver, const int matIdx, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const int idx);
   };

   template <typename TStepper> void SparseLinearCoordinatorImpl::buildSolverMatrixWrapper(SharedPtrMacro<TStepper > spSolver, const int matIdx, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const int idx)
   {
      // Resize operator matrix
      spSolver->rLHSMatrix(matIdx).resize(spEq->couplingInfo(comp).systemN(idx), spEq->couplingInfo(comp).systemN(idx));

      // Get model operator
      DecoupledZSparse  linOp;
      spEq->buildModelMatrix(linOp, ModelOperator::IMPLICIT_LINEAR, comp, idx, true);
      Solver::internal::addOperators(spSolver->rLHSMatrix(matIdx), 1.0, linOp); 
      
      // Solver is initialized
      spSolver->setInitialized();
   }
}
}

#endif // SPARSELINEARCOORDINATORIMPL_HPP
