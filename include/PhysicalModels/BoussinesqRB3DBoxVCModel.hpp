/** 
 * @file BoussinesqRB3DBoxVCModel.hpp
 * @brief Implementation of the Boussinesq Rayleigh-Benard 3D box (velocity-continuity formulation) model
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef BOUSSINESQRB3DBOXVCMODEL_HPP
#define BOUSSINESQRB3DBOXVCMODEL_HPP

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
#include "SpatialSchemes/3D/TTTScheme.hpp"

namespace GeoMHDiSCC {

   /**
    * @brief Implementation of the Boussinesq Rayleigh-Benard 3D box (velocity-continuity formulation) model
    */
   class BoussinesqRB3DBoxVCModel
   {
      public:
         /// Typedef for the spatial scheme used
         static const int DIMENSION = 3;

         /// Python script/module name
         static const std::string PYMODULE;

         /// Python model class name
         static const std::string PYCLASS;

         /// Typedef for the spatial scheme used
         typedef Schemes::TTTScheme SchemeType;

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
         BoussinesqRB3DBoxVCModel();

         /**
          * @brief Destructor
          */
         ~BoussinesqRB3DBoxVCModel();
   };

}

// 
// Block compilation of unusable parallelisation algorithms
//
#ifdef GEOMHDISCC_MPIALGO_SINGLE1D
#error "The SINGLE1D parallelisation is not supported!" 
#endif //GEOMHDISCC_MPIALGO_SINGLE1D
#ifdef GEOMHDISCC_MPIALGO_SINGLE2D
#error "The SINGLE2D parallelisation is not supported!" 
#endif //GEOMHDISCC_MPIALGO_SINGLE2D
#ifdef GEOMHDISCC_MPIALGO_TUBULAR
#error "The TUBULAR parallelisation is not supported!" 
#endif //GEOMHDISCC_MPIALGO_TUBULAR

#endif // BOUSSINESQRB3DBOXVCMODEL_HPP