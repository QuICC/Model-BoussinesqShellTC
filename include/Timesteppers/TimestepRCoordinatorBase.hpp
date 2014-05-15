/** 
 * @file TimestepRCoordinatorBase.hpp
 * @brief Implementation of a real fieldtimestep coordinator
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef TIMESTEPRCOORDINATORBASE_HPP
#define TIMESTEPRCOORDINATORBASE_HPP

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
    * @brief Implementation of real field timestepper
    */
   template <> class TimestepCoordinatorBase<false>: public Solver::SparseLinearCoordinatorBase<SparseTimestepper,false>
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
          * @brief Build the real operator, real field solver matrix
          *
          * @param spSolver   Shared sparse real solver
          * @param matIdx     Index of the solver matrix
          * @param spEq       Shared pointer to equation
          * @param comp       Field component
          * @param idx        Matrix index
          */
         virtual void buildSolverMatrix(TimestepCoordinatorBase::SharedRRSolverType spSolver, const int matIdx, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const int idx);

         /**
          * @brief Update time dependence
          */
         void updateTimeMatrices(const MHDFloat lhsCoeff, const MHDFloat rhsCoeff, const int step);

         /**
          * @brief Update the solvers
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

#endif // TIMESTEPRCOORDINATORBASE_HPP
