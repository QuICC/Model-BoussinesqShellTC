/** 
 * @file BoussinesqDynamoShellStdModel.hpp
 * @brief Implementation of the Boussinesq thermal convection dynamo in a spherical shell model (Toroidal/Poloidal formulation) without coupled solve (standard implementation)
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef BOUSSINESQDYNAMOSHELLSTDMODEL_HPP
#define BOUSSINESQDYNAMOSHELLSTDMODEL_HPP

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
// GEOMHDISCC_SPATIALSCHEME_FORMULATION = TORPOL;

namespace GeoMHDiSCC {

   /**
    * @brief Implementation of the Boussinesq thermal convection dynamo in a spherical shell model (Toroidal/Poloidal formulation) without coupled solve (standard implementation)
    */
   class BoussinesqDynamoShellStdModel
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
         BoussinesqDynamoShellStdModel();

         /**
          * @brief Destructor
          */
         ~BoussinesqDynamoShellStdModel();
   };

}

#endif // BOUSSINESQDYNAMOSHELLSTDMODEL_HPP
