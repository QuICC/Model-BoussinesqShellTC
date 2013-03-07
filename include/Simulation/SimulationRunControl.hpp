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
          */
         void update();

         /**
          * @brief Should the simulation keep running?
          */
         Runtime::Status::Id status() const;

         /**
          * @brief Update the status from control file
          */
         void checkFile();
         
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

      private:
   };
}

#endif // SIMULATIONRUNCONTROL_HPP
