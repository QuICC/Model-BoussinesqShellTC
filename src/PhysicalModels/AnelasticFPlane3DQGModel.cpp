/** \file AnelasticFPlane3DQGModel.cpp
 *  \brief Source of the anelastic f-plane 3DQG physical model
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "PhysicalModels/AnelasticFPlane3DQGModel.hpp"

// Project includes
//
#include "Enums/PhysicalNames.hpp"
#include "IoVariable/StateFileReader.hpp"
#include "IoVariable/StateFileWriter.hpp"
#include "IoVariable/VisualizationFileWriter.hpp"
#include "IoTools/IdToHuman.hpp"
#include "Equations/Asymptotics/AnelasticFPlane3DQG/AnelasticFPlane3DQGStreamfunction.hpp"
#include "Equations/Asymptotics/AnelasticFPlane3DQG/AnelasticFPlane3DQGVertical.hpp"
#include "Equations/Asymptotics/AnelasticFPlane3DQG/AnelasticFPlane3DQGTransport.hpp"

namespace GeoMHDiSCC {

   std::vector<PhysicalNames::Id> AnelasticFPlane3DQGModel::fieldIds()
   {
      // Create storage
      std::vector<PhysicalNames::Id> ids;

      // Add streamfunction
      ids.push_back(PhysicalNames::STREAMFUNCTION);

      // Add axial velocity
      ids.push_back(PhysicalNames::VELOCITYZ);

      // Add temperature fluctuations
      ids.push_back(PhysicalNames::TEMPERATURE);

      // Add mean temperature
      ids.push_back(PhysicalNames::MEANTEMPERATURE);

      return ids;
   }

   std::vector<std::string> AnelasticFPlane3DQGModel::boundaryNames()
   {
      // Create storage
      std::vector<std::string>   names;

      // Field IDs iterator
      std::vector<GeoMHDiSCC::PhysicalNames::Id>::const_iterator  it;
      std::vector<GeoMHDiSCC::PhysicalNames::Id> ids = AnelasticFPlane3DQGModel::fieldIds();

      for(it = ids.begin(); it != ids.end(); ++it)
      {
         names.push_back(IoTools::IdToHuman::toTag(*it));
      }

      return names;
   }

   std::vector<bool> AnelasticFPlane3DQGModel::isPeriodicBox()
   {
      std::vector<bool> box;

      // X direction is not periodic box
      box.push_back(false);

      // Y direction is periodic box
      box.push_back(true);

      // Z direction is not periodic box
      box.push_back(true);

      return box;
   }

   void AnelasticFPlane3DQGModel::addEquations(SharedSimulation spSim)
   {
      // Add mean temperature equation
      spSim->addScalarEquation<Equations::AnelasticFPlane3DQGTransport>();
      
      // Add transport equation
      spSim->addScalarEquation<Equations::AnelasticFPlane3DQGTransport>();
      
      // Add streamfunction equation
      spSim->addScalarEquation<Equations::AnelasticFPlane3DQGStreamfunction>();
      
      // Add vertical velocity equation
      spSim->addScalarEquation<Equations::AnelasticFPlane3DQGVertical>();
   }

   void AnelasticFPlane3DQGModel::addAsciiOutputFiles(SharedSimulation spSim)
   {
      // Add ASCII output file
      //pSim->addOutputFile(AN_ASCIIFILE);
      
      // Add ASCII output file
      //pSim->addOutputFile(AN_ASCIIFILE);
   }

   void AnelasticFPlane3DQGModel::addHdf5OutputFiles(SharedSimulation spSim)
   {
      // Field IDs iterator
      std::vector<GeoMHDiSCC::PhysicalNames::Id>::const_iterator  it;
      std::vector<GeoMHDiSCC::PhysicalNames::Id> ids = AnelasticFPlane3DQGModel::fieldIds();

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

   SharedSimulationBoundary AnelasticFPlane3DQGModel::createBoundary(const std::map<std::string,int>& bcIds)
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

   void AnelasticFPlane3DQGModel::setInitialState(SharedSimulation spSim)
   {
      // Field IDs iterator
      std::vector<GeoMHDiSCC::PhysicalNames::Id>::const_iterator  it;
      std::vector<GeoMHDiSCC::PhysicalNames::Id> ids = AnelasticFPlane3DQGModel::fieldIds();

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
