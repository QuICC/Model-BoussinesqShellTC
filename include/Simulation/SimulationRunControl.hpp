/** 
 * @file SimulationRunControl.hpp
 * @brief Implementation of a general simulation control structure
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef SIMULATIONRUNCONTROL_HPP
#define SIMULATIONRUNCONTROL_HPP

// Configuration includes
//
#include "SmartPointers/SharedPtrMacro.h"

// System includes
//

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Enums/RuntimeStatus.hpp"
#include "IoControl/ControlInterface.hpp"

namespace GeoMHDiSCC {

   /**
    * @brief Implementation of simulation control structure
    */
   class SimulationRunControl
   {
      public:
         /**
          * @brief Constructor
          */
         SimulationRunControl();

         /**
          * @brief Destructor
          */
         virtual ~SimulationRunControl();

         /**
          * @brief Update control status with simulation information
          *
          * @param simTime Simulation time
          * @param simDt   Simulation timestep
          */
         void updateSimulation(const MHDFloat simTime, const MHDFloat simDt);

         /**
          * @brief Update control status with cluster information
          *
          * @param wallTime Wall time
          */
         void updateCluster(const MHDFloat wallTime);

         /**
          * @brief Should the simulation keep running?
          */
         RuntimeStatus::Id status() const;

         /**
          * @brief Update the status from control file
          */
         void checkFile();

         /***
          * @brief Set the maximum simulation time
          *
          * @param maxTime New maximum simulation time
          */
         void setMaxSimTime(const MHDFloat maxTime);

         /***
          * @brief Set the maximum wall time
          *
          * @param maxTime New maximum wall time
          */
         void setMaxWallTime(const MHDFloat maxTime);
         
         /**
          * @brief Print run information
          *
          * @param stream  Output stream
          */
         void printInfo(std::ostream& stream);

         /**
          * @brief Static method to close simulation cleanly through signal
          */
         static void handleSignal(int signum);
         
      protected:
         /**
          * @brief Status to be switched by signal
          */
         static RuntimeStatus::Id   SIGNAL_STATUS;

         /**
          * @brief Current runtime status
          */
         RuntimeStatus::Id mStatus;

         /**
          * @brief External control file
          */
         IoControl::ControlInterface  mCtrlFile;

         /**
          * @brief Number of timesteps done
          */
         int   mSteps;

         /**
          * @brief Maximum simulation time
          */
         MHDFloat mMaxSimTime;

         /**
          * @brief Maximum wall time
          */
         MHDFloat mMaxWallTime;

      private:
         /**
          * @brief Initialise the signal handler
          */
         void initSignalHandler();
   };
}

#endif // SIMULATIONRUNCONTROL_HPP
