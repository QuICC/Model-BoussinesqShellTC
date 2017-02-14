/** 
 * @file ModelFactory.hpp
 * @brief Implementation of the physical model factory
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef MODELFACTORY_HPP
#define MODELFACTORY_HPP

// First include
//
#include "Python/PythonHeader.hpp"

// Configuration includes
//
#include "Framework/FrameworkMacro.h"

// System includes
//

// External includes
//

// Project includes
//
#include "Python/PythonModelWrapper.hpp"
#include "Simulation/Simulation.hpp"
#include "IoTools/IdToHuman.hpp"
#include "Model/PhysicalModelBase.hpp"

namespace QuICC {

   /**
    * @brief Implementation of the physical model factory 
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
      SharedSimulation  spSim(new Simulation);

      // Initialize the python model wrapper
      PythonModelWrapper::init();
      PythonModelWrapper::import(TModel::PYMODULE);
      PythonModelWrapper::createModel(TModel::PYCLASS);

      // Create list of field ID strings for boundary conditions
      std::vector<PhysicalNames::Id> fields = Model::PhysicalModelBase::fieldIds();
      std::vector<PhysicalNames::Id>::iterator fIt;
      std::vector<std::string>   bcNames;
      for(fIt = fields.begin(); fIt != fields.end(); ++fIt)
      {
         bcNames.push_back(IoTools::IdToHuman::toTag(*fIt));
      }

      // Create list of nondimensional ID strings for physical parameters
      std::vector<NonDimensional::Id> params = Model::PhysicalModelBase::paramIds();
      std::vector<NonDimensional::Id>::iterator pIt;
      std::vector<std::string>   ndNames;
      for(pIt = params.begin(); pIt != params.end(); ++pIt)
      {
         ndNames.push_back(IoTools::IdToHuman::toTag(*pIt));
      }

      // Add configuration file and parameters
      spSim->setConfiguration<TModel::DIMENSION>(TModel::SchemeType::type(), Model::PhysicalModelBase::isPeriodicBox(), bcNames, ndNames);

      // Initialise simulation
      spSim->initBase();

      // Initialise resolution
      spSim->initResolution<typename TModel::SchemeType>();

      // Add equations
      TModel::addEquations(spSim);

      // Add ASCII output files
      TModel::addAsciiOutputFiles(spSim);

      // Add HDF5 output files
      TModel::addHdf5OutputFiles(spSim);

      // Add statistics output files
      TModel::addStatsOutputFiles(spSim);

      // Set the boundary conditions
      SharedSimulationBoundary spBcs = spSim->createBoundary<TModel>();

      // Initialise the simulation
      spSim->init(spBcs);

      // Set initial state
      TModel::setInitialState(spSim);

      return spSim;
   }

}

#endif // MODELFACTORY_HPP
