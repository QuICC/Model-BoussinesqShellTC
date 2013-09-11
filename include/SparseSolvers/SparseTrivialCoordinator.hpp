/** 
 * @file SparseTrivialCoordinator.hpp
 * @brief Implementation of the base for a general sparse trivial solver coordinator
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 * @version 0.9.0
 * @date 2013-09-11
 */

#ifndef SPARSETRIVIALCOORDINATOR_HPP
#define SPARSETRIVIALCOORDINATOR_HPP

// Configuration includes
//
#include "SmartPointers/SharedPtrMacro.h"

// System includes
//

// External includes
//

// Project includes
//
#include "SparseSolvers/SparseCoordinatorBase.hpp"
#include "SparseSolvers/SparseDTrivialSolver.hpp"
#include "SparseSolvers/SparseZTrivialSolver.hpp"
#include "Equations/IScalarEquation.hpp"
#include "Equations/IVectorEquation.hpp"

namespace GeoMHDiSCC {

namespace Solver {

   /**
    * @brief Implementation of the base for a general sparse trivial solver coordinatore
    */
   class SparseTrivialCoordinator: public SparseCoordinatorBase<SharedSparseZTrivialSolver, SharedSparseDTrivialSolver>
   {
      public:
         /**
          * @brief Constructor
          */
         SparseTrivialCoordinator();

         /**
          * @brief Destructor
          */
         virtual ~SparseTrivialCoordinator();

         /**
          * @brief Initialise solver coordinator
          *
          * @param scalEq  Shared scalar equations
          * @param vectEq  Shared vector equations
          */
         virtual void init(const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq);

         /**
          * @brief Solve the equations
          *
          * @param scalEq  Shared scalar equations
          * @param vectEq  Shared vector equations
          * @param scalVar Shared scalar variables
          * @param vectVar Shared vector variables
          */
         void solve(const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq, const ScalarVariable_map& scalVar, const VectorVariable_map& vectVar);
         
      protected:
         /**
          * @brief Create a real linear solver
          */
         virtual void addSolverD(const int start);

         /**
          * @brief Create a complex linear solver
          */
         virtual void addSolverZ(const int start);

      private:
   };
}
}

#endif // SPARSETRIVIALCOORDINATORBASE_HPP
