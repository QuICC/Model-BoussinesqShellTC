/** \file SimulationRunControl.hpp
 *  \brief Implementation of a general simulation control structure
 *
 *  \mhdBug Needs test
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
#include "Enums/Runtime.hpp"
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
          * @brief Update control status
          *
          * @param simTime Simulation time
          * @param simDt   Simulation timestep
          */
         void update(const MHDFloat simTime, const MHDFloat simDt);

         /**
          * @brief Should the simulation keep running?
          */
         Runtime::Status::Id status() const;

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
         
      protected:
         /**
          * @brief Current runtime status
          */
         Runtime::Status::Id mStatus;

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
   };
}

#endif // SIMULATIONRUNCONTROL_HPP
