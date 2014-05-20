/** 
 * @file BoussinesqBeta3DQGModel.hpp
 * @brief Implementation of the Boussinesq beta 3DQG model
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef BOUSSINESQBETA3DQGMODEL_HPP
#define BOUSSINESQBETA3DQGMODEL_HPP

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
#include "SpatialSchemes/3D/TFTScheme.hpp"

namespace GeoMHDiSCC {

   /**
    * @brief Implementation of the Boussinesq beta 3DQG model
    */
   class BoussinesqBeta3DQGModel
   {
      public:
         /// Typedef for the spatial scheme used
         static const int DIMENSION = 3;

         /// Python script name
         static const std::string PYNAME;

         /// Typedef for the spatial scheme used
         typedef Schemes::TFTScheme SchemeType;

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
         BoussinesqBeta3DQGModel();

         /**
          * @brief Destructor
          */
         ~BoussinesqBeta3DQGModel();
   };

}

#endif // BOUSSINESQBETA3DQGMODEL_HPP
