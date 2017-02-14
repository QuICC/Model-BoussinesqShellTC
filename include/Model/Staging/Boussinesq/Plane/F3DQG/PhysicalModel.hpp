/** 
 * @file PhysicalModel.hpp
 * @brief Implementation of the Boussinesq F-plane 3DQG model
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef QUICC_MODEL_BOUSSINESQ_PLANE_F3DQG_PHYSICALMODEL_HPP
#define QUICC_MODEL_BOUSSINESQ_PLANE_F3DQG_PHYSICALMODEL_HPP

// Configuration includes
//

// System includes
//
#include <string>

// External includes
//

// Project includes
//
#include "Simulation/Simulation.hpp"
#include "Generator/StateGenerator.hpp"
#include "Generator/VisualizationGenerator.hpp"
#include "SpatialSchemes/3D/TFFScheme.hpp"

namespace QuICC {

namespace Model {

namespace Boussinesq {

namespace Plane {

namespace F3DQG {

   /**
    * @brief Implementation of the Boussinesq F-plane 3DQG model
    */
   class PhysicalModel
   {
      public:
         /// Typedef for the spatial scheme used
         static const int DIMENSION = 3;

         /// Python script/module name
         static const std::string PYMODULE;

         /// Python model class name
         static const std::string PYCLASS;

         /// Typedef for the spatial scheme used
         typedef Schemes::TFFScheme SchemeType;

         /**
          * @brief Add the required equations
          *
          * @param spSim   Shared simulation object
          */
         static void addEquations(SharedSimulation spSim);

         /**
          * @brief Add the initial state generation equations
          *
          * @param spGen   Shared generator object
          */
         static void addStates(SharedStateGenerator spGen);

         /**
          * @brief Add the visualization generation equations
          *
          * @param spGen   Shared visualization generator
          */
         static void addVisualizers(SharedVisualizationGenerator spVis);

         /**
          * @brief Set the visualization initial state
          *
          * @param spSim   Shared visualization generator
          */
         static void setVisualizationState(SharedVisualizationGenerator spVis);

         /**
          * @brief Add the required statistics output files
          *
          * @param spSim   Shared simulation object
          */
         static void addStatsOutputFiles(SharedSimulation spSim);

         /**
          * @brief Add the required ASCII output files
          *
          * @param spSim   Shared simulation object
          */
         static void addAsciiOutputFiles(SharedSimulation spSim);

         /**
          * @brief Add the required HDF5 output files
          *
          * @param spSim   Shared simulation object
          */
         static void addHdf5OutputFiles(SharedSimulation spSim);

         /**
          * @brief Set the initial state
          *
          * @param spSim   Shared simulation object
          */
         static void setInitialState(SharedSimulation spSim);

      protected:

      private:
         /**
          * @brief Constructor
          */
         PhysicalModel();

         /**
          * @brief Destructor
          */
         ~PhysicalModel();
   };

// 
// Block compilation of unusable timestepping schemes
//
#ifdef QUICC_TIMESTEPPER_IMEXRKCB2
#error "The ImExRKCB2 timestepper is not supported!" 
#endif //QUICC_TIMESTEPPER_IMEXRKCB2
#ifdef QUICC_TIMESTEPPER_IMEXRKCB3A
#error "The ImExRKCB3A timestepper is not supported!" 
#endif //QUICC_TIMESTEPPER_IMEXRKCB3A
#ifdef QUICC_TIMESTEPPER_IMEXRKCB3B
#error "The ImExRKCB3B timestepper is not supported!" 
#endif //QUICC_TIMESTEPPER_IMEXRKCB3B
#ifdef QUICC_TIMESTEPPER_IMEXRKCB3C
#error "The ImExRKCB3C timestepper is not supported!" 
#endif //QUICC_TIMESTEPPER_IMEXRKCB3C
#ifdef QUICC_TIMESTEPPER_IMEXRKCB3D
#error "The ImExRKCB3D timestepper is not supported!" 
#endif //QUICC_TIMESTEPPER_IMEXRKCB3D
#ifdef QUICC_TIMESTEPPER_IMEXRKCB3E
#error "The ImExRKCB3E timestepper is not supported!" 
#endif //QUICC_TIMESTEPPER_IMEXRKCB3E
#ifdef QUICC_TIMESTEPPER_IMEXRKCB3F
#error "The ImExRKCB3F timestepper is not supported!" 
#endif //QUICC_TIMESTEPPER_IMEXRKCB3F
#ifdef QUICC_TIMESTEPPER_IMEXRKCB4
#error "The ImExRKCB4 timestepper is not supported!" 
#endif //QUICC_TIMESTEPPER_IMEXRKCB4

}
}
}
}
}

#endif // QUICC_MODEL_BOUSSINESQ_PLANE_F3DQG_PHYSICALMODEL_HPP
