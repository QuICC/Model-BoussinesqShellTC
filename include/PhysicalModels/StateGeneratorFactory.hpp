/** 
 * @file StateGeneratorFactory.hpp
 * @brief Implementation of the state generator model factory
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef STATEGENERATORFACTORY_HPP
#define STATEGENERATORFACTORY_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "Generator/StateGenerator.hpp"
#include "IoTools/IdToHuman.hpp"

namespace GeoMHDiSCC {

   /**
    * @brief Implementation of the state generator model factory 
    */
   template <class TModel> class StateGeneratorFactory
   {
      public:
         /**
          * @brief Create a shared state generator for the model
          */
         static SharedStateGenerator createGenerator();

      protected:

      private:
         /**
          * @brief Constructor
          */
         StateGeneratorFactory();

         /**
          * @brief Destructor
          */
         ~StateGeneratorFactory();
   };

   template <class TModel> SharedStateGenerator StateGeneratorFactory<TModel>::createGenerator()
   {
      // Create simulation
      SharedStateGenerator  spGen(new StateGenerator);

      // Create list of field ID strings for boundary conditions
      std::vector<PhysicalNames::Id> fields = TModel::fieldIds();
      std::vector<PhysicalNames::Id>::iterator fIt;
      std::vector<std::string>   bcNames;
      for(fIt = fields.begin(); fIt != fields.end(); ++fIt)
      {
         bcNames.push_back(IoTools::IdToHuman::toTag(*fIt));
      }

      // Create list of nondimensional ID strings for physical parameters
      std::vector<NonDimensional::Id> params = TModel::paramIds();
      std::vector<NonDimensional::Id>::iterator pIt;
      std::vector<std::string>   ndNames;
      for(pIt = params.begin(); pIt != params.end(); ++pIt)
      {
         ndNames.push_back(IoTools::IdToHuman::toTag(*pIt));
      }

      // Add configuration file and parameters
      spGen->setConfiguration<TModel::DIMENSION>(TModel::SchemeType::type(), TModel::isPeriodicBox(), bcNames, ndNames);

      // Initialise simulation
      spGen->initBase();

      // Initialise resolution
      spGen->initResolution<typename TModel::SchemeType>();

      // Add initial state generation equations
      TModel::addStates(spGen);

      // Set the boundary conditions
      SharedSimulationBoundary spBcs = spGen->createBoundary<TModel>();

      // Initialise the simulation
      spGen->init(spBcs);

      return spGen;
   }

}

#endif // STATEGENERATORFACTORY_HPP
