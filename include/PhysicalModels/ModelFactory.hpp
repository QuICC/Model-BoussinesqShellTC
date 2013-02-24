/** \file ModelFactory.hpp
 *  \brief Implementation of an example physical model
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

      // Add equations
      TModel::addEquations(spSim);

      // Add configuration file
      TModel::setConfigurationFile(spSim);

      // Add initial state file
      TModel::setInitialStateFile(spSim);

      // Add ASCII output files
      TModel::addAsciiOutputFiles(spSim);

      // Add HDF5 output files
      TModel::addHdf5OutputFiles(spSim);

      return spSim;
   }

}

#endif // MODELFACTORY_HPP
