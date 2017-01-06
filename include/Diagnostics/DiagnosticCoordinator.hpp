/**
 * @file DiagnosticCoordinator.hpp
 * @brief Coordinator for the diagnostics computations 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef DIAGNOSTICCOORDINATOR_HPP
#define DIAGNOSTICCOORDINATOR_HPP

// Configuration includes
//

// System includes
//
#include <string>
#include <vector>

// External includes
//

// Project includes
//
#include "Enums/FieldIds.hpp"
#include "Diagnostics/ICflWrapper.hpp"
#include "TypeSelectors/VariableSelector.hpp"

namespace QuICC {

namespace Diagnostics {

   /**
    * @brief Coordinator for the diagnostics computations
    */
   class DiagnosticCoordinator
   {
      public:
         /**
          * @brief Constructor
          */
         DiagnosticCoordinator();

         /**
          * @brief Constructor
          */
         ~DiagnosticCoordinator();

         /**
          * @brief Initialise the coordinator
          *
          * @param mesh    Vector of grid values
          * @param scalars Map of shared scalar variables
          * @param vectors Map of shared vector variables
          * @param tstep   Timestep information
          */
         void init(const std::vector<Array>& mesh, const std::map<PhysicalNames::Id, Datatypes::SharedScalarVariableType>&  scalars, const std::map<PhysicalNames::Id, Datatypes::SharedVectorVariableType>&  vectors, const Array& tstep); 

         /**
          * @brief Compute the initial CFL condition
          */
         void initialCfl();

         /**
          * @brief Compute the current CFL condition
          */
         void updateCfl();

         /**
          * @brief Get CFL condition
          */
         MHDFloat cfl() const;

         /**
          * @brief Get max error goal
          */
         MHDFloat maxError() const;

         /**
          * @brief Get start time
          */
         MHDFloat startTime() const;

         /**
          * @brief Get start timestep
          */
         MHDFloat startTimestep() const;

         /**
          * @brief Use time and timestep from state file
          */
         void useStateTime(const MHDFloat time, const MHDFloat timestep);

         /**
          * @brief Synchronize diagnostics among CPUs
          */
         void synchronize();

      protected:

      private:
         /**
          * @brief Maximum timestep
          */
         const MHDFloat mcMaxStep;

         /**
          * @brief Minimum timestep
          */
         const MHDFloat mcMinStep;

         /**
          * @brief Fixed timestep
          */
         MHDFloat mFixedStep;

         /**
          * @brief Max error goal
          */
         MHDFloat mMaxError;

         /**
          * @brief Current CFL condition
          */
         MHDFloat mCfl;

         /**
          * @brief Start simulation time
          */
         MHDFloat mStartTime;

         /**
          * @brief Start simulation timestep
          */
         MHDFloat mStartTimestep;

         /**
          * @brief Shared pointer to a CFL condition wrapper
          */
         SharedICflWrapper  mspCflWrapper;
   };
}
}

#endif // DIAGNOSTICCOORDINATOR_HPP
