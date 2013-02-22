/** \file SimulationBase.hpp
 *  \brief Building block for the implementation of a simulation
 */

#ifndef SIMULATIONBASE_HPP
#define SIMULATIONBASE_HPP

// Configuration includes
//
#include "Simulation/PrepMacros/EquationParametersMacro.h"

// System includes
//
#include <map>

// External includes
//

// Project includes
//
#include "Simulation/Enums/PhysicalNames.hpp"
#include "Simulation/System/SpatialSimulationBase.hpp"
#include "Simulation/IO/IOSystem.hpp"
#include "Simulation/Controls/SimulationControl.hpp"
#include "Simulation/Timesteppers/Timestepper.hpp"

namespace EPMPhoenix {

   /**
    * \brief Building block for the implementation of a simulation
    */
   class SimulationBase: public SpatialSimulationBase
   {
      public:
         /**
          * @brief Simple empty destructor
          */
         virtual ~SimulationBase() {};

      protected:
         /**
          * @brief Constructor
          */
         SimulationBase();

         /**
          * @brief Initialise the simulation base
          */
         virtual void initSimulation();

         /**
          * @brief Initialise the equations
          */
         virtual void initEquations() = 0;

         /**
          * @brief Initialise the input files
          */
         virtual void initOutput() = 0;

         /**
          * @brief Init base ouput files
          */
         void initBaseOutput();

         /**
          * @brief Setup the ouput files
          */
         void setupOutput();

         /**
          * @brief Initialise the timestepper
          */
         void initTimestepper();

         /**
          * @brief Initialise the input files
          */
         void initInput();

         /**
          * @brief Cleanup unused memory from simulation base
          */
         virtual void cleanupSimulation();

         /**
          * @brief Do operations required just before starting the time integration
          */
         void preRun();

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

         /**
          * @brief Finalise the output files
          */
         void finaliseOutput();

         /**
          * @brief Finalise the simulation base
          */
         void finaliseSimulation();

         /**
          * @brief Get the equation parameters
          */
         SharedEquationParameters spEqParams();

         /**
          * @brief Get/Set the IO system
          */
         IOSystem& ioSys();

         /**
          * @brief Get/Set the simulation control
          */
         SimulationControl& simCtrl();

         /**
          * @brief Get/Set the timestepper
          */
         Timestepper& timestepper();

      private:
         /**
          * @brief IO system
          */
         IOSystem mIOSystem;

         /**
          * @brief Shared simulation control
          */
         SimulationControl   mSimControl;

         /**
          * @brief Storage for the shared equation parameters
          */
         SharedEquationParameters   mspEqParams;

         /**
          * @brief Timestepper 
          */
         Timestepper mTimestepper;
   };

   inline SharedEquationParameters SimulationBase::spEqParams()
   {
      return this->mspEqParams;
   }

   inline IOSystem& SimulationBase::ioSys()
   {
      return this->mIOSystem;
   }

   inline SimulationControl& SimulationBase::simCtrl()
   {
      return this->mSimControl;
   }

   inline Timestepper& SimulationBase::timestepper()
   {
      return this->mTimestepper;
   }
}

#endif // SIMULATIONBASE_HPP
