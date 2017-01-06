/** 
 * @file BoussinesqRTCShellModel.hpp
 * @brief Implementation of the Boussinesq rotating thermal convection in a spherical shell (Toroidal/Poloidal formulation)
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef BOUSSINESQRTCSHELLMODEL_HPP
#define BOUSSINESQRTCSHELLMODEL_HPP

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
#include "SpatialSchemes/3D/SLFmScheme.hpp"

// THIS IS NOT A COMMENT BUT AND OPTION READ BY CMAKE
// QUICC_SPATIALSCHEME_FORMULATION = TORPOL;

namespace GeoMHDiSCC {

   /**
    * @brief Implementation of the Boussinesq rotating thermal convection spherical shell model (Toroidal/Poloidal formulation)
    */
   class BoussinesqRTCShellModel
   {
      public:
         /// Typedef for the spatial scheme used
         static const int DIMENSION = 3;

         /// Python script/module name
         static const std::string PYMODULE;

         /// Python model class name
         static const std::string PYCLASS;

         /// Typedef for the spatial scheme used
         typedef Schemes::SLFmScheme SchemeType;

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
         BoussinesqRTCShellModel();

         /**
          * @brief Destructor
          */
         ~BoussinesqRTCShellModel();
   };

}

// 
// Block compilation of unusable parallelisation algorithms
//
#ifdef QUICC_MPIALGO_SINGLE1D
#error "The SINGLE1D parallelisation is not supported!" 
#endif //QUICC_MPIALGO_SINGLE1D
#if defined QUICC_MPIALGO_TUBULAR && !defined QUICC_SPLINALG_MUMPS && !defined QUICC_MPISPSOLVE
#error "The TUBULAR parallelisation is not supported!" 
#endif //defined QUICC_MPIALGO_TUBULAR && !defined QUICC_SPLINALG_MUMPS && !defined QUICC_MPISPSOLVE

#endif // BOUSSINESQRTCSHELLMODEL_HPP
