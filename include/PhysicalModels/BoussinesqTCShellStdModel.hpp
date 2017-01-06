/** 
 * @file BoussinesqTCShellStdModel.hpp
 * @brief Implementation of the Boussinesq thermal convection in a spherical shell (Toroidal/Poloidal formulation) without coupled solve (standard implementation)
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef BOUSSINESQTCSHELLSTDMODEL_HPP
#define BOUSSINESQTCSHELLSTDMODEL_HPP

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
#include "SpatialSchemes/3D/SLFlScheme.hpp"

// THIS IS NOT A COMMENT BUT AND OPTION READ BY CMAKE
// QUICC_SPATIALSCHEME_FORMULATION = TORPOL;

namespace QuICC {

   /**
    * @brief Implementation of the Boussinesq thermal convection spherical shell model (Toroidal/Poloidal formulation) without coupled solve (standard implementation)
    */
   class BoussinesqTCShellStdModel
   {
      public:
         /// Typedef for the spatial scheme used
         static const int DIMENSION = 3;

         /// Python script/module name
         static const std::string PYMODULE;

         /// Python model class name
         static const std::string PYCLASS;

         /// Typedef for the spatial scheme used
         typedef Schemes::SLFlScheme SchemeType;

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
          * @brief Add the required statistics output files
          * 
          * @param spSim   Shared simulation object
          */
         static void addStatsOutputFiles(SharedSimulation spSim);

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
         BoussinesqTCShellStdModel();

         /**
          * @brief Destructor
          */
         ~BoussinesqTCShellStdModel();
   };

}

#endif // BOUSSINESQTCSHELLSTDMODEL_HPP
