/** \file BoussinesqBetaSphGModel.cpp
 *  \brief Source of the Boussinesq beta 3DQG model with spherical gravity
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "PhysicalModels/BoussinesqBetaSphGModel.hpp"

// Project includes
//
#include "Enums/FieldIds.hpp"
#include "IoVariable/StateFileReader.hpp"
#include "IoVariable/StateFileWriter.hpp"
#include "IoVariable/VisualizationFileWriter.hpp"
#include "IoTools/IdToHuman.hpp"
#include "Equations/Asymptotics/Beta3DQG/BoussinesqBetaSphGStreamfunction.hpp"
#include "Equations/Asymptotics/Beta3DQG/BoussinesqBetaSphGVertical.hpp"
#include "Equations/Asymptotics/Beta3DQG/BoussinesqBetaSphGTransport.hpp"

namespace GeoMHDiSCC {

   std::vector<PhysicalNames::Id> BoussinesqBetaSphGModel::fieldIds()
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

   std::vector<std::string> BoussinesqBetaSphGModel::boundaryNames()
   {
      // Create storage
      std::vector<std::string>   names;

      // Field IDs iterator
      std::vector<GeoMHDiSCC::PhysicalNames::Id>::const_iterator  it;
      std::vector<GeoMHDiSCC::PhysicalNames::Id> ids = BoussinesqBetaSphGModel::fieldIds();

      for(it = ids.begin(); it != ids.end(); ++it)
      {
         names.push_back(IoTools::IdToHuman::toTag(*it));
      }

      return names;
   }

   std::vector<bool> BoussinesqBetaSphGModel::isPeriodicBox()
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

   void BoussinesqBetaSphGModel::addEquations(SharedSimulation spSim)
   {
      // Add transport equation
      spSim->addScalarEquation<Equations::BoussinesqBetaSphGTransport>();
      
      // Add streamfunction equation
      spSim->addScalarEquation<Equations::BoussinesqBetaSphGStreamfunction>();
      
      // Add vertical velocity equation
      spSim->addScalarEquation<Equations::BoussinesqBetaSphGVertical>();
   }

   void BoussinesqBetaSphGModel::addAsciiOutputFiles(SharedSimulation spSim)
   {
      // Add ASCII output file
      //pSim->addOutputFile(AN_ASCIIFILE);
      
      // Add ASCII output file
      //pSim->addOutputFile(AN_ASCIIFILE);
   }

   void BoussinesqBetaSphGModel::addHdf5OutputFiles(SharedSimulation spSim)
   {
      // Field IDs iterator
      std::vector<GeoMHDiSCC::PhysicalNames::Id>::const_iterator  it;
      std::vector<GeoMHDiSCC::PhysicalNames::Id> ids = BoussinesqBetaSphGModel::fieldIds();

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

   SharedSimulationBoundary BoussinesqBetaSphGModel::createBoundary(const std::map<std::string,int>& bcIds)
   {
      // Create shared simulation boundary
      SharedSimulationBoundary  spBcs(new SimulationBoundary());

      // Storage for the dimension ID
      Dimensions::Simulation::Id dimId;

      // Create equation and field keys
      SpectralFieldId eqId;
      SpectralFieldId fieldId;

      // Temperature equation
      //    ... boundary conditions
      eqId = std::make_pair(PhysicalNames::TEMPERATURE, FieldComponents::Spectral::SCALAR);
      spBcs->initStorage(eqId);
      dimId = Dimensions::Simulation::SIM1D;
      fieldId = std::make_pair(PhysicalNames::TEMPERATURE, FieldComponents::Spectral::SCALAR);
      spBcs->initBcStorage(eqId, fieldId, dimId);
      if(bcIds.find(IoTools::IdToHuman::toTag(eqId.first))->second == 0)
      {
         spBcs->addBc(eqId, fieldId, dimId, Spectral::BoundaryConditions::VALUE, Spectral::IBoundary::LEFT);
         spBcs->addBc(eqId, fieldId, dimId, Spectral::BoundaryConditions::VALUE, Spectral::IBoundary::RIGHT);
      } else
      {
         throw Exception("Unknown temperature boundary conditions in configuration file");
      }

      // Streamfunction equation
      //    ... boundary conditions
      eqId = std::make_pair(PhysicalNames::STREAMFUNCTION, FieldComponents::Spectral::SCALAR);
      spBcs->initStorage(eqId);
      dimId = Dimensions::Simulation::SIM1D;
      fieldId = std::make_pair(PhysicalNames::STREAMFUNCTION, FieldComponents::Spectral::SCALAR);
      spBcs->initBcStorage(eqId, fieldId, dimId);
      spBcs->addBc(eqId, fieldId, dimId, Spectral::BoundaryConditions::VALUE, Spectral::IBoundary::LEFT); 
      spBcs->addBc(eqId, fieldId, dimId, Spectral::BoundaryConditions::VALUE, Spectral::IBoundary::RIGHT); 
      // No-slip boundary conditions
      if(bcIds.find(IoTools::IdToHuman::toTag(eqId.first))->second == 0)
      {
         spBcs->addBc(eqId, fieldId, dimId, Spectral::BoundaryConditions::FIRST_DERIVATIVE, Spectral::IBoundary::LEFT); 
         spBcs->addBc(eqId, fieldId, dimId, Spectral::BoundaryConditions::FIRST_DERIVATIVE, Spectral::IBoundary::RIGHT); 
      // Stress-free boundary conditions
      } else if(bcIds.find(IoTools::IdToHuman::toTag(eqId.first))->second == 1)
      {
         spBcs->addBc(eqId, fieldId, dimId, Spectral::BoundaryConditions::SECOND_DERIVATIVE, Spectral::IBoundary::LEFT); 
         spBcs->addBc(eqId, fieldId, dimId, Spectral::BoundaryConditions::SECOND_DERIVATIVE, Spectral::IBoundary::RIGHT);
      } else
      {
         throw Exception("Unknown streamfunction boundary conditions in configuration file");
      }
      dimId = Dimensions::Simulation::SIM3D;
      spBcs->initBcStorage(eqId, fieldId, dimId);
      spBcs->addBc(eqId, fieldId, dimId, Spectral::BoundaryConditions::BETA_SLOPE, Spectral::IBoundary::LEFT); 
      //    ... coupled boundary conditions
      fieldId = std::make_pair(PhysicalNames::VELOCITYZ, FieldComponents::Spectral::SCALAR);
      dimId = Dimensions::Simulation::SIM3D;
      spBcs->initBcStorage(eqId, fieldId, dimId); 
      spBcs->addBc(eqId, fieldId, dimId, Spectral::BoundaryConditions::VALUE, Spectral::IBoundary::LEFT); 

      // Axial velocity equation
      //    ... boundary conditions
      eqId = std::make_pair(PhysicalNames::VELOCITYZ, FieldComponents::Spectral::SCALAR);
      spBcs->initStorage(eqId);
      dimId = Dimensions::Simulation::SIM1D;
      fieldId = std::make_pair(PhysicalNames::VELOCITYZ, FieldComponents::Spectral::SCALAR);
      spBcs->initBcStorage(eqId, fieldId, dimId);
      // No-slip boundary conditions
      if(bcIds.find(IoTools::IdToHuman::toTag(eqId.first))->second == 0)
      {
         spBcs->addBc(eqId, fieldId, dimId, Spectral::BoundaryConditions::VALUE, Spectral::IBoundary::LEFT); 
         spBcs->addBc(eqId, fieldId, dimId, Spectral::BoundaryConditions::VALUE, Spectral::IBoundary::RIGHT); 
      // Stress-free boundary conditions
      } else if(bcIds.find(IoTools::IdToHuman::toTag(eqId.first))->second == 1)
      {
         spBcs->addBc(eqId, fieldId, dimId, Spectral::BoundaryConditions::FIRST_DERIVATIVE, Spectral::IBoundary::LEFT); 
         spBcs->addBc(eqId, fieldId, dimId, Spectral::BoundaryConditions::FIRST_DERIVATIVE, Spectral::IBoundary::RIGHT); 
      } else
      {
         throw Exception("Unknown velocityz boundary conditions in configuration file");
      }
      dimId = Dimensions::Simulation::SIM3D;
      spBcs->initBcStorage(eqId, fieldId, dimId);
      spBcs->addBc(eqId, fieldId, dimId, Spectral::BoundaryConditions::VALUE, Spectral::IBoundary::RIGHT); 
      //    ... coupled boundary conditions
      fieldId = std::make_pair(PhysicalNames::STREAMFUNCTION, FieldComponents::Spectral::SCALAR);
      dimId = Dimensions::Simulation::SIM3D;
      spBcs->initBcStorage(eqId, fieldId, dimId); 
      spBcs->addBc(eqId, fieldId, dimId, Spectral::BoundaryConditions::BETA_SLOPE, Spectral::IBoundary::RIGHT); 

      return spBcs;
   }

   void BoussinesqBetaSphGModel::setInitialState(SharedSimulation spSim)
   {
      // Field IDs iterator
      std::vector<GeoMHDiSCC::PhysicalNames::Id>::const_iterator  it;
      std::vector<GeoMHDiSCC::PhysicalNames::Id> ids = BoussinesqBetaSphGModel::fieldIds();

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
