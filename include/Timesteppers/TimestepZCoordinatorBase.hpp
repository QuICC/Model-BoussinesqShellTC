/** 
 * @file TimestepZCoordinatorBase.hpp
 * @brief Implementation of a complex field timestep coordinator
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef TIMESTEPZCOORDINATORBASE_HPP
#define TIMESTEPZCOORDINATORBASE_HPP

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
#include "SparseSolvers/SparseLinearCoordinatorBase.hpp"
#include "Timesteppers/SparseTimestepper.hpp"
#include "Equations/IScalarEquation.hpp"
#include "Equations/IVectorEquation.hpp"

namespace GeoMHDiSCC {

namespace Timestep {

   // Forward declaration
   template <bool TIsComplex> class TimestepCoordinatorBase;

   /**
    * @brief Implementation of real field timestep coordinator
    */
   template <> class TimestepCoordinatorBase<true>: public Solver::SparseLinearCoordinatorBase<SparseTimestepper,true>
   {
      public:
         /**
          * @brief Constructor
          */
         TimestepCoordinatorBase();

         /**
          * @brief Destructor
          */
         ~TimestepCoordinatorBase();
         
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
         virtual void buildSolverMatrix(TimestepCoordinatorBase::SharedRZSolverType spSolver, const int matIdx, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const int idx);

         /**
          * @brief Build the complex operator, complex field solver matrix
          *
          * @param spSolver   Shared sparse complex solver
          * @param matIdx     Index of the solver matrix
          * @param spEq       Shared pointer to equation
          * @param comp       Field component
          * @param idx        Matrix index
          */
         virtual void buildSolverMatrix(TimestepCoordinatorBase::SharedZZSolverType spSolver, const int matIdx, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const int idx);

         /**
          * @brief Update time dependence
          */
         void updateTimeMatrices(const MHDFloat lhsCoeff, const MHDFloat rhsCoeff, const int step);

         /**
          * @brief Update solvers
          */
         void updateSolvers();

         /**
          * @brief Compute the RHS of all linear systems
          */
         void computeRHS();

      private:
   };
}
}

#endif // TIMESTEPZCOORDINATORBASE_HPP
