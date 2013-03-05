/** \file ModelFactory.hpp
 *  \brief Implementation of an example physical model
 *
 *  \mhdBug Needs test
 */

#ifndef MODELFACTORY_HPP
#define MODELFACTORY_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "Simulation/Simulation.hpp"

namespace GeoMHDiSCC {

   /**
    * @brief Implementation of an example physical model
    */
   template <class TModel> class ModelFactory
   {
      public:
         /**
          * @brief Create a shared simulation for the model
          */
         static SharedSimulation createSimulation();

      protected:

      private:
         /**
          * @brief Constructor
          */
         ModelFactory();

         /**
          * @brief Destructor
          */
         ~ModelFactory();
   };

   template <class TModel> SharedSimulation ModelFactory<TModel>::createSimulation()
   {
      // Create simulation
      SharedSimulation  spSim;

      // Add configuration file and parameters
      spSim->setConfiguration<TModel::DIMENSION, typename TModel::ParametersType>(TModel::boundaryNames());

      // Initialise simulation
      spSim->initBase();

      // Add equations
      spSim->initResolution<typename TModel::SchemeType>();

      // Add equations
      TModel::addEquations(spSim);

      // Add ASCII output files
      TModel::addAsciiOutputFiles(spSim);

      // Add HDF5 output files
      TModel::addHdf5OutputFiles(spSim);

      // Set the boundary conditions
      SharedPtrMacro<SimulationBoundary> spBcs = spSim->createBoundary<TModel>();

      // Initialise the simulation
      spSim->init(*spBcs);

      // Set initial state
      TModel::setInitialState(spSim);

      return spSim;
   }

}

#endif // MODELFACTORY_HPP
