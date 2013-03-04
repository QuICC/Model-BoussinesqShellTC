/** \file ExampleModel.hpp
 *  \brief Implementation of an example physical model
 */

#ifndef EXAMPLEMODEL_HPP
#define EXAMPLEMODEL_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "Simulation/Simulation.hpp"
#include "SpatialSchemes/3D/TFTScheme.hpp"
#include "Equations/Parameters/PrRaXGParameters.hpp"

namespace GeoMHDiSCC {

   /**
    * @brief Implementation of an example physical model
    *
    *  \mhdBug Needs implementation
    *  \mhdBug Needs test
    */
   class ExampleModel
   {
      public:
         /// Typedef for the spatial scheme used
         static const int DIMENSION = 3;

         /// Typedef for the spatial scheme used
         typedef TFTScheme SchemeType;

         /// Typedef for the equation parameters used
         typedef Equations::PrRaXGParameters  ParametersType;

         /**
          * @brief Get the required equations
          *
          * @param spSim   Shared simulation object
          */
         static void addEquations(SharedSimulation spSim);

         /**
          * @brief Get the required XML input files
          *
          * @param spSim   Shared simulation object
          */
         static void setConfigurationFile(SharedSimulation spSim);

         /**
          * @brief Set the initial state
          *
          * @param spSim   Shared simulation object
          */
         static void setInitialState(SharedSimulation spSim);

         /**
          * @brief Get the required ASCII output files
          *
          * @param spSim   Shared simulation object
          */
         static void addAsciiOutputFiles(SharedSimulation spSim);

         /**
          * @brief Get the required HDF5 output files
          *
          * @param spSim   Shared simulation object
          */
         static void addHdf5OutputFiles(SharedSimulation spSim);

      protected:

      private:
         /**
          * @brief Constructor
          */
         ExampleModel();

         /**
          * @brief Simple empty destructor
          */
         ~ExampleModel();
   };

}

#endif // EXAMPLEMODEL_HPP
