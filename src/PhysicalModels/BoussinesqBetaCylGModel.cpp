/** \file BoussinesqBetaCylGModel.cpp
 *  \brief Source of the Boussinesq beta 3DQG model with cylindrical gravity
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "PhysicalModels/BoussinesqBetaCylGModel.hpp"

// Project includes
//
#include "Enums/FieldIds.hpp"
#include "IoVariable/StateFileReader.hpp"
#include "IoVariable/StateFileWriter.hpp"
#include "IoVariable/VisualizationFileWriter.hpp"
#include "IoTools/IdToHuman.hpp"
#include "Equations/Asymptotics/Beta3DQG/Boussinesq/BoussinesqBetaCylGStreamfunction.hpp"
#include "Equations/Asymptotics/Beta3DQG/Boussinesq/BoussinesqBetaCylGVertical.hpp"
#include "Equations/Asymptotics/Beta3DQG/Boussinesq/BoussinesqBetaCylGTransport.hpp"
#include "Equations/Asymptotics/Beta3DQG/Boussinesq/BoussinesqBetaCylGVorticity.hpp"
#include "Generator/States/RandomScalarState.hpp"
#include "Generator/States/ExactScalarState.hpp"
#include "Generator/Visualizers/FieldVisualizer.hpp"

namespace GeoMHDiSCC {

   std::vector<PhysicalNames::Id> BoussinesqBetaCylGModel::fieldIds()
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

   std::vector<NonDimensional::Id> BoussinesqBetaCylGModel::paramIds()
   {
      // Create storage
      std::vector<NonDimensional::Id> ids;

      // Add Prandtl number
      ids.push_back(NonDimensional::PRANDTL);

      // Add Rayleigh number
      ids.push_back(NonDimensional::RAYLEIGH);

      // Add gamma
      ids.push_back(NonDimensional::GAMMA);

      // Add chi
      ids.push_back(NonDimensional::CHI);

      return ids;
   }

   std::vector<bool> BoussinesqBetaCylGModel::isPeriodicBox()
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

   void BoussinesqBetaCylGModel::addEquations(SharedSimulation spSim)
   {
      // Add transport equation
      spSim->addScalarEquation<Equations::BoussinesqBetaCylGTransport>();
      
      // Add streamfunction equation
      spSim->addScalarEquation<Equations::BoussinesqBetaCylGStreamfunction>();
      
      // Add vertical velocity equation
      spSim->addScalarEquation<Equations::BoussinesqBetaCylGVertical>();
      
      // Add vorticity computation
      spSim->addScalarEquation<Equations::BoussinesqBetaCylGVorticity>();
   }

   void BoussinesqBetaCylGModel::addStates(SharedStateGenerator spGen)
   {
      // Shared pointer to equation
      Equations::SharedRandomScalarState spRand;
      // Shared pointer to equation
      Equations::SharedExactScalarState spExact;

      // Add transport initial state generation equation
      spRand = spGen->addScalarEquation<Equations::RandomScalarState>();
      spRand->setIdentity(PhysicalNames::TEMPERATURE);
      spRand->setSpectrum(-10,10, 1e4, 1e4, 1e4);
      
      // Add streamfunction initial state generation equation
      spRand = spGen->addScalarEquation<Equations::RandomScalarState>();
      spRand->setIdentity(PhysicalNames::STREAMFUNCTION);
      spRand->setSpectrum(-10,10, 1e4, 1e4, 1e4);
      
      // Add vertical velocity initial state generation equation
      spExact = spGen->addScalarEquation<Equations::ExactScalarState>();
      spExact->setIdentity(PhysicalNames::VELOCITYZ);
      spExact->setStateType(Equations::ExactScalarState::SINE);

      // Add output file
      IoVariable::SharedStateFileWriter spOut(new IoVariable::StateFileWriter(SchemeType::type(), SchemeType::isRegular()));
      spOut->expect(PhysicalNames::TEMPERATURE);
      spOut->expect(PhysicalNames::STREAMFUNCTION);
      spOut->expect(PhysicalNames::VELOCITYZ);
      spGen->addOutputFile(spOut);
   }

   void BoussinesqBetaCylGModel::addVisualizers(SharedVisualizationGenerator spVis)
   {
      // Shared pointer to basic field visualizer
      Equations::SharedFieldVisualizer spField;

      // Add transport field visualization
      spField = spVis->addScalarEquation<Equations::FieldVisualizer>();
      spField->setFields(true, false);
      spField->setIdentity(PhysicalNames::TEMPERATURE);
      
      // Add streamfunction field visualization
      spField = spVis->addScalarEquation<Equations::FieldVisualizer>();
      spField->setFields(true, false);
      spField->setIdentity(PhysicalNames::STREAMFUNCTION);
      
      // Add vertical velocity field visualization
      spField = spVis->addScalarEquation<Equations::FieldVisualizer>();
      spField->setFields(true, false);
      spField->setIdentity(PhysicalNames::VELOCITYZ);

      // Add output file
      IoVariable::SharedVisualizationFileWriter spOut(new IoVariable::VisualizationFileWriter(SchemeType::type()));
      spOut->expect(PhysicalNames::TEMPERATURE);
      spOut->expect(PhysicalNames::STREAMFUNCTION);
      spOut->expect(PhysicalNames::VELOCITYZ);
      spVis->addOutputFile(spOut);
   }

   void BoussinesqBetaCylGModel::setVisualizationState(SharedVisualizationGenerator spVis)
   {
      // Create and add initial state file to IO
      IoVariable::SharedStateFileReader spIn(new IoVariable::StateFileReader("4Visu", SchemeType::type(), SchemeType::isRegular()));

      // Set expected fields
      spIn->expect(PhysicalNames::TEMPERATURE);
      spIn->expect(PhysicalNames::STREAMFUNCTION);
      spIn->expect(PhysicalNames::VELOCITYZ);

      // Set simulation state
      spVis->setInitialState(spIn);
   }

   void BoussinesqBetaCylGModel::addAsciiOutputFiles(SharedSimulation spSim)
   {
      // Add ASCII output file
      //pSim->addOutputFile(AN_ASCIIFILE);
      
      // Add ASCII output file
      //pSim->addOutputFile(AN_ASCIIFILE);
   }

   void BoussinesqBetaCylGModel::addHdf5OutputFiles(SharedSimulation spSim)
   {
      // Field IDs iterator
      std::vector<GeoMHDiSCC::PhysicalNames::Id>::const_iterator  it;
      std::vector<GeoMHDiSCC::PhysicalNames::Id> ids = BoussinesqBetaCylGModel::fieldIds();

      // Create and add state file to IO
      IoVariable::SharedStateFileWriter spState(new IoVariable::StateFileWriter(SchemeType::type(), SchemeType::isRegular()));
      for(it = ids.begin(); it != ids.end(); ++it)
      {
         spState->expect(*it);
      }
      spSim->addOutputFile(spState);
   }

   SharedSimulationBoundary BoussinesqBetaCylGModel::createBoundary(const std::map<std::string,int>& bcIds)
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

   void BoussinesqBetaCylGModel::setInitialState(SharedSimulation spSim)
   {
      // Field IDs iterator
      std::vector<GeoMHDiSCC::PhysicalNames::Id>::const_iterator  it;
      std::vector<GeoMHDiSCC::PhysicalNames::Id> ids = BoussinesqBetaCylGModel::fieldIds();

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
