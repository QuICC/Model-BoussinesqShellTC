/** 
 * @file BoussinesqRBCAnnulusVCModel.hpp
 * @brief Implementation of Boussinesq Rayleigh-Benard convection in a annulus (velocity-continuity formulation)
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef BOUSSINESQRBCANNULUSVCMODEL_HPP
#define BOUSSINESQRBCANNULUSVCMODEL_HPP

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
#include "SpatialSchemes/3D/AFTScheme.hpp"

namespace QuICC {

   /**
    * @brief Implementation of Boussinesq Rayleigh-Benard convection in a annulus (velocity-continuity formulation)
    */
   class BoussinesqRBCAnnulusVCModel
   {
      public:
         /// Typedef for the spatial scheme used
         static const int DIMENSION = 3;

         /// Python script/module name
         static const std::string PYMODULE;

         /// Python model class name
         static const std::string PYCLASS;

         /// Typedef for the spatial scheme used
         typedef Schemes::AFTScheme SchemeType;

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
         BoussinesqRBCAnnulusVCModel();

         /**
          * @brief Destructor
          */
         ~BoussinesqRBCAnnulusVCModel();
   };

}

// 
// Block compilation of unusable parallelisation algorithms
//
#ifdef QUICC_MPIALGO_SINGLE1D
#error "The SINGLE1D parallelisation is not supported!" 
#endif //QUICC_MPIALGO_SINGLE1D
#ifdef QUICC_MPIALGO_SINGLE2D
#error "The SINGLE2D parallelisation is not supported!" 
#endif //QUICC_MPIALGO_SINGLE2D
#ifdef QUICC_MPIALGO_TUBULAR
#error "The TUBULAR parallelisation is not supported!" 
#endif //QUICC_MPIALGO_TUBULAR

#endif // BOUSSINESQRBCANNULUSVCMODEL_HPP