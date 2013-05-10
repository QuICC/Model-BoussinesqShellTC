/** \file RotatingRBModel.hpp
 *  \brief Implementation of the rotating Rayleigh-Benard physical model
 */

#ifndef ROTATINGRBMODEL_HPP
#define ROTATINGRBMODEL_HPP

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
#include "SpatialSchemes/3D/TFFScheme.hpp"

namespace GeoMHDiSCC {

   /**
    * @brief Implementation of the rotating Rayleigh-Benard physical model
    */
   class RotatingRBModel
   {
      public:
         /// Typedef for the spatial scheme used
         static const int DIMENSION = 3;

         /// Typedef for the spatial scheme used
         typedef TFFScheme SchemeType;

         /**
          * @brief Get vector of names for the boundary conditions
          */
         static std::vector<PhysicalNames::Id> fieldIds();

         /**
          * @brief Get vector of names for the nondimensional parameters
          */
         static std::vector<NonDimensional::Id> paramIds();

         /**
          * @brief Get vector of bools about periodic box
          */
         static std::vector<bool> isPeriodicBox();

         /**
          * @brief Add the required equations
          *
          * @param spSim   Shared simulation object
          */
         static void addEquations(SharedSimulation spSim);

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
          * @brief Set the boundary conditions
          *
          * @param bcIds Boundary condition IDs
          */
         static SharedSimulationBoundary createBoundary(const std::map<std::string, int>& bcIds);

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
         RotatingRBModel();

         /**
          * @brief Destructor
          */
         ~RotatingRBModel();
   };

}

#endif // ROTATINGRBMODEL_HPP
