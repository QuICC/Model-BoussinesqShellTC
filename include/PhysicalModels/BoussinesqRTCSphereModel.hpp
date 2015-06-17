/** 
 * @file BoussinesqRTCSphereModel.hpp
 * @brief Implementation of the Boussinesq rotating thermal convection in a sphere (Toroidal/Poloidal formulation)
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef BOUSSINESQRTCSPHEREMODEL_HPP
#define BOUSSINESQRTCSPHEREMODEL_HPP

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
#include "SpatialSchemes/3D/BLFmScheme.hpp"

// THIS IS NOT A COMMENT BUT AND OPTION READ BY CMAKE
// GEOMHDISCC_SPATIALSCHEME_FORMULATION = TORPOL;

namespace GeoMHDiSCC {

   /**
    * @brief Implementation of the Boussinesq rotating thermal convection sphere model (Toroidal/Poloidal formulation)
    */
   class BoussinesqRTCSphereModel
   {
      public:
         /// Typedef for the spatial scheme used
         static const int DIMENSION = 3;

         /// Python script/module name
         static const std::string PYMODULE;

         /// Python model class name
         static const std::string PYCLASS;

         /// Typedef for the spatial scheme used
         typedef Schemes::BLFmScheme SchemeType;

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
         BoussinesqRTCSphereModel();

         /**
          * @brief Destructor
          */
         ~BoussinesqRTCSphereModel();
   };

}

// 
// Block compilation of unusable parallelisation algorithms
//
#ifdef GEOMHDISCC_MPIALGO_SINGLE1D
#error "The SINGLE1D parallelisation is not supported!" 
#endif //GEOMHDISCC_MPIALGO_SINGLE1D
#if defined GEOMHDISCC_MPIALGO_TUBULAR && !defined GEOMHDISCC_SPLINALG_MUMPS && !defined GEOMHDISCC_MPISPSOLVE
#error "The TUBULAR parallelisation is not supported!" 
#endif //defined GEOMHDISCC_MPIALGO_TUBULAR && !defined GEOMHDISCC_SPLINALG_MUMPS && !defined GEOMHDISCC_MPISPSOLVE

#endif // BOUSSINESQRTCSPHEREMODEL_HPP