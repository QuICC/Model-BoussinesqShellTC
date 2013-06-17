/** \file Simulation.hpp
 *  \brief High level implementation of a simulation
 */

#ifndef SIMULATION_HPP
#define SIMULATION_HPP

// Configuration includes
//
#include "SmartPointers/SharedPtrMacro.h"

// System includes
//

// External includes
//

// Project includes
//
#include "Simulation/SimulationBase.hpp"
#include "Timesteppers/TimestepCoordinator.hpp"

namespace GeoMHDiSCC {

   /**
    * @brief High level implementation of a simulation's execution steps.
    */
   class Simulation: public SimulationBase
   {
      public:
         /**
          * @brief Constructor
          */
         Simulation();

         /**
          * @brief Simple empty destructor
          */
         virtual ~Simulation();

      protected:

      private:
         /**
          * @brief Initialise the generator specific base components
          */
         virtual void initAdditionalBase();

         /**
          * @brief Do operations required just before starting the time integration
          */
         virtual void preRun();

         /**
          * @brief Do operations required during the main loop
          */
         virtual void mainRun();

         /**
          * @brief Do operations required just after finishing the time integration
          */
         virtual void postRun();

         /**
          * @brief Solve all equations
          */
         void solveEquations();

         /**
          * @brief Timestep the prognostic equations
          */
         void solvePrognosticEquations();

         /**
          * @brief Write the output if required
          */
         void writeOutput();

         /**
          * @brief Timestep coordinator
          */
         Timestep::TimestepCoordinator mTimestepCoordinator;
   };

   /// Typedef for a shared pointer of a Simulation
   typedef SharedPtrMacro<Simulation> SharedSimulation;

}

#endif // SIMULATION_HPP
