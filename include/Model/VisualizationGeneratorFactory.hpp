/** 
 * @file VisualizationGeneratorFactory.hpp
 * @brief Implementation of the visualization generator model factory
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef VISUALIZATIONGENERATORFACTORY_HPP
#define VISUALIZATIONGENERATORFACTORY_HPP

// First include
//
#include "Python/PythonHeader.hpp"

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "Python/PythonModelWrapper.hpp"
#include "Generator/VisualizationGenerator.hpp"
#include "IoTools/IdToHuman.hpp"
#include "Model/PhysicalModelBase.hpp"

namespace QuICC {

   /**
    * @brief Implementation of the visualization generator model factory 
    */
   template <class TModel> class VisualizationGeneratorFactory
   {
      public:
         /**
          * @brief Create a shared state generator for the model
          */
         static SharedVisualizationGenerator createVisualization();

      protected:

      private:
         /**
          * @brief Constructor
          */
         VisualizationGeneratorFactory();

         /**
          * @brief Destructor
          */
         ~VisualizationGeneratorFactory();
   };

   template <class TModel> SharedVisualizationGenerator VisualizationGeneratorFactory<TModel>::createVisualization()
   {
      // Create simulation
      SharedVisualizationGenerator  spVis(new VisualizationGenerator);

      // Initialize the python model wrapper
      PythonModelWrapper::init();
      PythonModelWrapper::import(TModel::PYMODULE);

      // Try to setup specialized Visualization model
      PythonModelWrapper::createModel(TModel::PYCLASS, "Visu");

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
      spVis->setConfiguration<TModel::DIMENSION>(TModel::SchemeType::type(), Model::PhysicalModelBase::isPeriodicBox(), bcNames, ndNames);

      // Initialise simulation
      spVis->initBase();

      // Initialise resolution
      spVis->initResolution<typename TModel::SchemeType>();

      // Add initial state generation equations
      TModel::addVisualizers(spVis);

      // Set the boundary conditions
      SharedSimulationBoundary spBcs = spVis->createBoundary<TModel>();

      // Initialise the simulation
      spVis->init(spBcs);

      // Set initial state
      TModel::setVisualizationState(spVis);

      return spVis;
   }

}

#endif // VISUALIZATIONGENERATORFACTORY_HPP
