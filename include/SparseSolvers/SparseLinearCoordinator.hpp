/** \file SparseLinearCoordinator.hpp
 *  \brief Implementation of a general sparse linear solver coordinator
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
#include "SparseSolvers/SparseLinearCoordinatorBase.hpp"
#include "SparseSolvers/SparseDLinearSolver.hpp"
#include "SparseSolvers/SparseZLinearSolver.hpp"
#include "Equations/IScalarEquation.hpp"
#include "Equations/IVectorEquation.hpp"

namespace GeoMHDiSCC {

namespace Solver {

   /**
    * @brief Implementation of general sparse linear solver coordinatore
    */
   class SparseLinearCoordinator: public SparseLinearCoordinatorBase
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
          */
         void solve(const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq);
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
          * @param isLhs      Flag to update LHS and RHS time dependent matrix
          */
         virtual void buildSolverMatrix(SharedSparseDLinearSolver spSolver, const int matIdx, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const int idx, const bool isLhs);

         /**
          * @brief Build the complex solver matrix
          *
          * @param spSolver   Shared sparse real solver
          * @param matIdx     Index of the solver matrix
          * @param spEq       Shared pointer to equation
          * @param comp       Field component
          * @param idx        Matrix index
          * @param isLhs      Flag to update LHS and RHS time dependent matrix
          */
         virtual void buildSolverMatrix(SharedSparseDLinearSolver spSolver, const int matIdx, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const int idx, const bool isLhs);

      private:
   };
}
}

#endif // SPARSELINEARCOORDINATOR_HPP
