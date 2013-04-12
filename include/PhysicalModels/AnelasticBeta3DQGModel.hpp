/** \file AnelasticBeta3DQGModel.hpp
 *  \brief Implementation of the anelastic beta 3DQG physical model
 */

#ifndef ANELASTICBETA3DQGMODEL_HPP
#define ANELASTICBETA3DQGMODEL_HPP

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
#include "SpatialSchemes/3D/TFTScheme.hpp"
#include "Equations/Parameters/PrRaXGParameters.hpp"

namespace GeoMHDiSCC {

   /**
    * @brief Implementation of the anelastic beta 3DQG physical model
    */
   class AnelasticBeta3DQGModel
   {
      public:
         /// Typedef for the spatial scheme used
         static const int DIMENSION = 3;

         /// Typedef for the spatial scheme used
         typedef TFTScheme SchemeType;

         /// Typedef for the equation parameters used
         typedef Equations::PrRaXGParameters  ParametersType;

         /**
          * @brief Get vector of names for the boundary conditions
          */
         static std::vector<PhysicalNames::Id> fieldIds();

         /**
          * @brief Get vector of names for the boundary conditions
          */
         static std::vector<std::string> boundaryNames();

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
         static SharedPtrMacro<SimulationBoundary> createBoundary(const std::map<std::string, int>& bcIds);

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
         AnelasticBeta3DQGModel();

         /**
          * @brief Destructor
          */
         ~AnelasticBeta3DQGModel();
   };

}

#endif // ANELASTICBETA3DQGMODEL_HPP