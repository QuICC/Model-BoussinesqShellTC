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
#include "Timers/ExecutionTimer.hpp"
#include "Simulation/SimulationRunControl.hpp"
#include "Simulation/SimulationIoControl.hpp"
#include "Timesteppers/Timestepper.hpp"
#include "IoConfig/ConfigurationReader.hpp"

namespace GeoMHDiSCC {

   /**
    * @brief High level implementation of a simulation's execution steps.
    */
   class Simulation
   {
      public:
         /**
          * @brief Constructor
          *
          * The constructor simply calls the constructor of TSimImpl.
          */
         Simulation();

         /**
          * @brief Simple empty destructor
          */
         ~Simulation();

         /**
          * @brief Initialise the different components of the simulation
          */
         void init();

         /**
          * @brief Run the simulation
          */
         void run();

         /**
          * @brief Finalise simulation run
          */
         void finalize();

         /**
          * @brief Add scalar equation to solver
          *
          * @param spEq Shared scalar equation
          *
          * \mhdBug Fake implementation
          */
         void addEquation(int spEq);//SharedScalarEquation  spEq);

         /**
          * @brief Add vector equation to solver
          *
          * @param spEq Shared vector equation
          *
          * \mhdBug Fake implementation
          */
         void addEquation(double spEq);//SharedVectorEquation  spEq);

         /**
          * @brief Set the simulation configuration file
          *
          * @param spCfgFile Shared configuration file
          *
          * \mhdBug Fake implementation
          */
         void setConfigurationFile(SharedConfigurationFile spCfgFile);

         /**
          * @brief Set solver initial state file
          *
          * @param spInitFile Shared initial state file
          *
          * \mhdBug Fake implementation
          */
         void setInitialStateFile(int spInitFile);//Shared spInitFile);

         /**
          * @brief Add ASCII output file to solver
          *
          * @param spOutFile Shared ASCII output file
          *
          * \mhdBug Fake implementation
          */
         void addOutputFile(int spOutFile);//SharedAscii spOutFile);

         /**
          * @brief Add HDF5 output file to solver
          *
          * @param spOutFile Shared HDF5 output file
          *
          * \mhdBug Fake implementation
          */
         void addOutputFile(double spOutFile);//SharedHdf5 spOutFile);

      protected:
         /**
          * @brief Do operations required just before starting the time integration
          */
         void preRun();

         /**
          * @brief Compute the nonlinear terms
          */
         void computeNonlinear();

         /**
          * @brief Timestep the equations
          */
         void timestepEquations();

         /**
          * @brief Write the output if required
          */
         void writeOutput();

         /**
          * @brief Do operations required just after finishing the time integration
          */
         void postRun();

      private:
         /**
          * @brief Execution timer
          */
         ExecutionTimer mExecutionTimer;

         /**
          * @brief Simulation run control
          */
         SimulationRunControl   mSimRunCtrl;

         /**
          * @brief Timestepper
          */
         Timestepper mTimestepper;

         /**
          * @brief Simulation IO control
          */
         SimulationIoControl mSimIoCtrl;
   };

   /// Typedef for a shared pointer of a Simulation
   typedef SharedPtrMacro<Simulation> SharedSimulation;

}

#endif // SIMULATION_HPP
