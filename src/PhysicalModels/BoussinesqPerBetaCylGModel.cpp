/** 
 * @file BoussinesqPerBetaCylGModel.cpp
 * @brief Source of the Boussinesq beta 3DQG model with cylindrical gravity
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 * @version 0.9.0
 * @date 2013-09-11
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "PhysicalModels/BoussinesqPerBetaCylGModel.hpp"

// Project includes
//
#include "Enums/FieldIds.hpp"
#include "IoVariable/StateFileReader.hpp"
#include "IoVariable/StateFileWriter.hpp"
#include "IoVariable/VisualizationFileWriter.hpp"
#include "IoTools/IdToHuman.hpp"
#include "Equations/Asymptotics/Beta3DQG/Boussinesq/BoussinesqPerBetaCylGStreamfunction.hpp"
#include "Equations/Asymptotics/Beta3DQG/Boussinesq/BoussinesqPerBetaCylGVertical.hpp"
#include "Equations/Asymptotics/Beta3DQG/Boussinesq/BoussinesqPerBetaCylGTransport.hpp"

namespace GeoMHDiSCC {

   std::vector<PhysicalNames::Id> BoussinesqPerBetaCylGModel::fieldIds()
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

   std::vector<NonDimensional::Id> BoussinesqPerBetaCylGModel::paramIds()
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

   std::vector<bool> BoussinesqPerBetaCylGModel::isPeriodicBox()
   {
      std::vector<bool> box;

      // X direction is periodic box
      box.push_back(true);

      // Y direction is periodic box
      box.push_back(true);

      // Z direction is not periodic box
      box.push_back(false);

      return box;
   }

   void BoussinesqPerBetaCylGModel::addEquations(SharedSimulation spSim)
   {
      // Add transport equation
      spSim->addScalarEquation<Equations::BoussinesqPerBetaCylGTransport>();
      
      // Add streamfunction equation
      spSim->addScalarEquation<Equations::BoussinesqPerBetaCylGStreamfunction>();
      
      // Add vertical velocity equation
      spSim->addScalarEquation<Equations::BoussinesqPerBetaCylGVertical>();
   }

   void BoussinesqPerBetaPerCylGModel::addStates(SharedStateGenerator spGen)
   {
      throw Exception("Not implemented yet!");
   }

   void BoussinesqPerBetaCylGModel::addVisualizers(SharedVisualizationGenerator spVis)
   {
      throw Exception("Not implemented yet!");
   }

   void BoussinesqPerBetaCylGModel::setVisualizationState(SharedVisualizationGenerator spVis)
   {
      throw Exception("Not implemented yet!");
   }

   void BoussinesqPerBetaCylGModel::addAsciiOutputFiles(SharedSimulation spSim)
   {
      // Add ASCII output file
      //pSim->addOutputFile(AN_ASCIIFILE);
      
      // Add ASCII output file
      //pSim->addOutputFile(AN_ASCIIFILE);
   }

   void BoussinesqPerBetaCylGModel::addHdf5OutputFiles(SharedSimulation spSim)
   {
      // Field IDs iterator
      std::vector<GeoMHDiSCC::PhysicalNames::Id>::const_iterator  it;
      std::vector<GeoMHDiSCC::PhysicalNames::Id> ids = BoussinesqPerBetaCylGModel::fieldIds();

      // Create and add state file to IO
      IoVariable::SharedStateFileWriter spState(new IoVariable::StateFileWriter(SchemeType::type(), SchemeType::isRegular()));
      for(it = ids.begin(); it != ids.end(); ++it)
      {
         spState->expect(*it);
      }
      spSim->addOutputFile(spState);
   }

   SharedSimulationBoundary BoussinesqPerBetaCylGModel::createBoundary(const std::map<std::string,int>& bcIds)
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
      // No boundary condition in periodic radius case
      if(bcIds.find(IoTools::IdToHuman::toTag(eqId.first))->second != -1)
      {
         throw Exception("Unknown temperature boundary conditions in configuration file");
      }

      // Streamfunction equation
      //    ... boundary conditions
      eqId = std::make_pair(PhysicalNames::STREAMFUNCTION, FieldComponents::Spectral::SCALAR);
      if(bcIds.find(IoTools::IdToHuman::toTag(eqId.first))->second != -1)
      {
         throw Exception("Unknown streamfunction boundary conditions in configuration file");
      }
      spBcs->initStorage(eqId);
      dimId = Dimensions::Simulation::SIM1D;
      fieldId = std::make_pair(PhysicalNames::STREAMFUNCTION, FieldComponents::Spectral::SCALAR);
      spBcs->initBcStorage(eqId, fieldId, dimId);
      spBcs->addBc(eqId, fieldId, dimId, Spectral::BoundaryConditions::BETA_SLOPE, Spectral::IBoundary::LEFT); 
      //    ... coupled boundary conditions
      fieldId = std::make_pair(PhysicalNames::VELOCITYZ, FieldComponents::Spectral::SCALAR);
      dimId = Dimensions::Simulation::SIM1D;
      spBcs->initBcStorage(eqId, fieldId, dimId); 
      spBcs->addBc(eqId, fieldId, dimId, Spectral::BoundaryConditions::VALUE, Spectral::IBoundary::LEFT); 

      // Axial velocity equation
      //    ... boundary conditions
      eqId = std::make_pair(PhysicalNames::VELOCITYZ, FieldComponents::Spectral::SCALAR);
      // No-slip boundary conditions
      if(bcIds.find(IoTools::IdToHuman::toTag(eqId.first))->second != -1)
      {
         throw Exception("Unknown velocityz boundary conditions in configuration file");
      }
      spBcs->initStorage(eqId);
      dimId = Dimensions::Simulation::SIM1D;
      fieldId = std::make_pair(PhysicalNames::VELOCITYZ, FieldComponents::Spectral::SCALAR);
      spBcs->initBcStorage(eqId, fieldId, dimId);
      spBcs->addBc(eqId, fieldId, dimId, Spectral::BoundaryConditions::VALUE, Spectral::IBoundary::RIGHT); 
      //    ... coupled boundary conditions
      fieldId = std::make_pair(PhysicalNames::STREAMFUNCTION, FieldComponents::Spectral::SCALAR);
      dimId = Dimensions::Simulation::SIM1D;
      spBcs->initBcStorage(eqId, fieldId, dimId); 
      spBcs->addBc(eqId, fieldId, dimId, Spectral::BoundaryConditions::BETA_SLOPE, Spectral::IBoundary::RIGHT); 

      return spBcs;
   }

   void BoussinesqPerBetaCylGModel::setInitialState(SharedSimulation spSim)
   {
      // Field IDs iterator
      std::vector<GeoMHDiSCC::PhysicalNames::Id>::const_iterator  it;
      std::vector<GeoMHDiSCC::PhysicalNames::Id> ids = BoussinesqPerBetaCylGModel::fieldIds();

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
