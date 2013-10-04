/** 
 * @file SparseTrivialCoordinator.hpp
 * @brief Implementation of the base for a general sparse trivial solver coordinator
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
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
#include "SparseSolvers/SparseTrivialSolver.hpp"
#include "Equations/IScalarEquation.hpp"
#include "Equations/IVectorEquation.hpp"

namespace GeoMHDiSCC {

namespace Solver {

   /**
    * @brief Implementation of the base for a general sparse trivial solver coordinator
    */
   class SparseTrivialCoordinator: public SparseCoordinatorBase<SparseZTrivialSolver, SparseRZTrivialSolver>
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

      private:
   };
}
}

#endif // SPARSETRIVIALCOORDINATORBASE_HPP
