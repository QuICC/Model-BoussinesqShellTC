/** \file AnelasticBetaSphGModel.cpp
 *  \brief Source of the anelastic beta 3DQG model with spherical gravity
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "PhysicalModels/AnelasticBetaSphGModel.hpp"

// Project includes
//
#include "Enums/FieldIds.hpp"
#include "IoVariable/StateFileReader.hpp"
#include "IoVariable/StateFileWriter.hpp"
#include "IoVariable/VisualizationFileWriter.hpp"
#include "IoTools/IdToHuman.hpp"
#include "Equations/Asymptotics/AnelasticBeta3DQG/AnelasticBeta3DQGStreamfunction.hpp"
#include "Equations/Asymptotics/AnelasticBeta3DQG/AnelasticBeta3DQGVertical.hpp"
#include "Equations/Asymptotics/AnelasticBeta3DQG/AnelasticBeta3DQGTransport.hpp"

namespace GeoMHDiSCC {

   std::vector<PhysicalNames::Id> AnelasticBetaSphGModel::fieldIds()
   {
      // Create storage
      std::vector<PhysicalNames::Id> ids;

      // Add streamfunction
      ids.push_back(PhysicalNames::STREAMFUNCTION);

      // Add axial velocity
      ids.push_back(PhysicalNames::VELOCITYZ);

      // Add temperature
      ids.push_back(PhysicalNames::TEMPERATURE);

      return ids;
   }

   std::vector<std::string> AnelasticBetaSphGModel::boundaryNames()
   {
      // Create storage
      std::vector<std::string>   names;

      // Field IDs iterator
      std::vector<GeoMHDiSCC::PhysicalNames::Id>::const_iterator  it;
      std::vector<GeoMHDiSCC::PhysicalNames::Id> ids = AnelasticBetaSphGModel::fieldIds();

      for(it = ids.begin(); it != ids.end(); ++it)
      {
         names.push_back(IoTools::IdToHuman::toTag(*it));
      }

      return names;
   }

   std::vector<bool> AnelasticBetaSphGModel::isPeriodicBox()
   {
      std::vector<bool> box;

      // X direction is not periodic box
      box.push_back(false);

      // Y direction is periodic box
      box.push_back(true);

      // Z direction is not periodic box
      box.push_back(false);

      return box;
   }

   void AnelasticBetaSphGModel::addEquations(SharedSimulation spSim)
   {
      // Add transport equation
      spSim->addScalarEquation<Equations::AnelasticBeta3DQGTransport>();
      
      // Add streamfunction equation
      spSim->addScalarEquation<Equations::AnelasticBeta3DQGStreamfunction>();
      
      // Add vertical velocity equation
      spSim->addScalarEquation<Equations::AnelasticBeta3DQGVertical>();
   }

   void AnelasticBetaSphGModel::addAsciiOutputFiles(SharedSimulation spSim)
   {
      // Add ASCII output file
      //pSim->addOutputFile(AN_ASCIIFILE);
      
      // Add ASCII output file
      //pSim->addOutputFile(AN_ASCIIFILE);
   }

   void AnelasticBetaSphGModel::addHdf5OutputFiles(SharedSimulation spSim)
   {
      // Field IDs iterator
      std::vector<GeoMHDiSCC::PhysicalNames::Id>::const_iterator  it;
      std::vector<GeoMHDiSCC::PhysicalNames::Id> ids = AnelasticBetaSphGModel::fieldIds();

      // Create and add state file to IO
      IoVariable::SharedStateFileWriter spState(new IoVariable::StateFileWriter(SchemeType::type(), SchemeType::isRegular()));
      for(it = ids.begin(); it != ids.end(); ++it)
      {
         spState->expect(*it);
      }
      spSim->addOutputFile(spState);
      
      // Create and add visualization file to IO
      IoVariable::SharedVisualizationFileWriter spViz(new IoVariable::VisualizationFileWriter(SchemeType::type()));
      for(it = ids.begin(); it != ids.end(); ++it)
      {
         spViz->expect(*it);
      }
      spSim->addOutputFile(spViz);
   }

   SharedSimulationBoundary AnelasticBetaSphGModel::createBoundary(const std::map<std::string,int>& bcIds)
   {
      // Create shared simulation boundary
      SharedSimulationBoundary  spBcs(new SimulationBoundary());

      // Create boundary condition key
      Equations::IEvolutionEquation::BcKeyType  key;

      // Temperature equation
      //    ... boundary conditions
      spBcs->initStorage(PhysicalNames::TEMPERATURE);
      key = std::make_pair(FieldComponents::Spectral::SCALAR, Dimensions::Simulation::SIM1D);
      spBcs->initBcStorage(PhysicalNames::TEMPERATURE, key);
      if(bcIds.find(IoTools::IdToHuman::toTag(PhysicalNames::TEMPERATURE))->second == 0)
      {
         spBcs->addBc(PhysicalNames::TEMPERATURE, key, Spectral::BoundaryConditions::VALUE, Spectral::IBoundary::LEFT);
         spBcs->addBc(PhysicalNames::TEMPERATURE, key, Spectral::BoundaryConditions::VALUE, Spectral::IBoundary::RIGHT);
      } else
      {
         throw Exception("Unknown temperature boundary conditions in configuration file");
      }

      // Streamfunction equation
      //    ... boundary conditions
      spBcs->initStorage(PhysicalNames::STREAMFUNCTION);
      key = std::make_pair(FieldComponents::Spectral::SCALAR, Dimensions::Simulation::SIM1D);
      spBcs->initBcStorage(PhysicalNames::STREAMFUNCTION, key);
      spBcs->addBc(PhysicalNames::STREAMFUNCTION, key, Spectral::BoundaryConditions::VALUE, Spectral::IBoundary::LEFT); 
      spBcs->addBc(PhysicalNames::STREAMFUNCTION, key, Spectral::BoundaryConditions::VALUE, Spectral::IBoundary::RIGHT); 
      // No-slip boundary conditions
      if(bcIds.find(IoTools::IdToHuman::toTag(PhysicalNames::STREAMFUNCTION))->second == 0)
      {
         spBcs->addBc(PhysicalNames::STREAMFUNCTION, key, Spectral::BoundaryConditions::FIRST_DERIVATIVE, Spectral::IBoundary::LEFT); 
         spBcs->addBc(PhysicalNames::STREAMFUNCTION, key, Spectral::BoundaryConditions::FIRST_DERIVATIVE, Spectral::IBoundary::RIGHT); 
      // Stress-free boundary conditions
      } else if(bcIds.find(IoTools::IdToHuman::toTag(PhysicalNames::STREAMFUNCTION))->second == 1)
      {
         spBcs->addBc(PhysicalNames::STREAMFUNCTION, key, Spectral::BoundaryConditions::SECOND_DERIVATIVE, Spectral::IBoundary::LEFT); 
         spBcs->addBc(PhysicalNames::STREAMFUNCTION, key, Spectral::BoundaryConditions::SECOND_DERIVATIVE, Spectral::IBoundary::RIGHT);
      } else
      {
         throw Exception("Unknown streamfunction boundary conditions in configuration file");
      }
      key = std::make_pair(FieldComponents::Spectral::SCALAR, Dimensions::Simulation::SIM3D);
      spBcs->initBcStorage(PhysicalNames::STREAMFUNCTION, key);
      spBcs->addBc(PhysicalNames::STREAMFUNCTION, key, Spectral::BoundaryConditions::BETA_SLOPE, Spectral::IBoundary::LEFT); 
      //    ... coupled boundary conditions
      key = std::make_pair(FieldComponents::Spectral::SCALAR, Dimensions::Simulation::SIM3D);
      spBcs->initCBcStorage(PhysicalNames::STREAMFUNCTION, PhysicalNames::VELOCITYZ, key); 
      spBcs->addCBc(PhysicalNames::STREAMFUNCTION, PhysicalNames::VELOCITYZ, key, Spectral::BoundaryConditions::VALUE, Spectral::IBoundary::LEFT); 

      // Axial velocity equation
      //    ... boundary conditions
      spBcs->initStorage(PhysicalNames::VELOCITYZ);
      key = std::make_pair(FieldComponents::Spectral::SCALAR, Dimensions::Simulation::SIM1D);
      spBcs->initBcStorage(PhysicalNames::VELOCITYZ, key);
      // No-slip boundary conditions
      if(bcIds.find(IoTools::IdToHuman::toTag(PhysicalNames::VELOCITYZ))->second == 0)
      {
         spBcs->addBc(PhysicalNames::VELOCITYZ, key, Spectral::BoundaryConditions::VALUE, Spectral::IBoundary::LEFT); 
         spBcs->addBc(PhysicalNames::VELOCITYZ, key, Spectral::BoundaryConditions::VALUE, Spectral::IBoundary::RIGHT); 
      // Stress-free boundary conditions
      } else if(bcIds.find(IoTools::IdToHuman::toTag(PhysicalNames::VELOCITYZ))->second == 1)
      {
         spBcs->addBc(PhysicalNames::VELOCITYZ, key, Spectral::BoundaryConditions::FIRST_DERIVATIVE, Spectral::IBoundary::LEFT); 
         spBcs->addBc(PhysicalNames::VELOCITYZ, key, Spectral::BoundaryConditions::FIRST_DERIVATIVE, Spectral::IBoundary::RIGHT); 
      } else
      {
         throw Exception("Unknown velocityz boundary conditions in configuration file");
      }
      key = std::make_pair(FieldComponents::Spectral::SCALAR, Dimensions::Simulation::SIM3D);
      spBcs->initBcStorage(PhysicalNames::VELOCITYZ, key);
      spBcs->addBc(PhysicalNames::VELOCITYZ, key, Spectral::BoundaryConditions::VALUE, Spectral::IBoundary::RIGHT); 
      //    ... coupled boundary conditions
      key = std::make_pair(FieldComponents::Spectral::SCALAR, Dimensions::Simulation::SIM3D);
      spBcs->initCBcStorage(PhysicalNames::VELOCITYZ, PhysicalNames::STREAMFUNCTION, key); 
      spBcs->addCBc(PhysicalNames::VELOCITYZ, PhysicalNames::STREAMFUNCTION, key, Spectral::BoundaryConditions::BETA_SLOPE, Spectral::IBoundary::RIGHT); 

      return spBcs;
   }

   void AnelasticBetaSphGModel::setInitialState(SharedSimulation spSim)
   {
      // Field IDs iterator
      std::vector<GeoMHDiSCC::PhysicalNames::Id>::const_iterator  it;
      std::vector<GeoMHDiSCC::PhysicalNames::Id> ids = AnelasticBetaSphGModel::fieldIds();

      // Create and add initial state file to IO
      IoVariable::SharedStateFileReader spInit(new IoVariable::StateFileReader("_initial", SchemeType::type(), SchemeType::isRegular()));

      // Set expected field names
      for(it = ids.begin(); it != ids.end(); ++it)
      {
         spInit->expect(*it);
      }

      // Set simulation state
      spSim->setInitialState(spInit);
   }

}
