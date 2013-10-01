/** 
 * @file TestTTTModel.cpp
 * @brief Source of the test model for the TTT scheme
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "PhysicalModels/TestTTTModel.hpp"

// Project includes
//
#include "Enums/FieldIds.hpp"
#include "IoVariable/StateFileReader.hpp"
#include "IoVariable/StateFileWriter.hpp"
#include "IoVariable/VisualizationFileWriter.hpp"
#include "IoTools/IdToHuman.hpp"
#include "Generator/States/RandomScalarState.hpp"
#include "Generator/States/ExactScalarState.hpp"
#include "Generator/Visualizers/FieldVisualizer.hpp"

namespace GeoMHDiSCC {

   std::vector<PhysicalNames::Id> TestTTTModel::fieldIds()
   {
      // Create storage
      std::vector<PhysicalNames::Id> ids;

      // Add first scalar
      ids.push_back(PhysicalNames::TEMPERATURE);

      // Add first scalar
      ids.push_back(PhysicalNames::STREAMFUNCTION);

      // Add first scalar
      ids.push_back(PhysicalNames::VORTICITY);

      return ids;
   }

   std::vector<NonDimensional::Id> TestTTTModel::paramIds()
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

   std::vector<bool> TestTTTModel::isPeriodicBox()
   {
      std::vector<bool> box;

      // X direction is not periodic box
      box.push_back(false);

      // Y direction is periodic box
      box.push_back(false);

      // Z direction is not periodic box
      box.push_back(false);

      return box;
   }

   void TestTTTModel::addEquations(SharedSimulation spSim)
   {
   }

   void TestTTTModel::addStates(SharedStateGenerator spGen)
   {
      // Shared pointer to equation
      Equations::SharedRandomScalarState spRand;
      Equations::SharedExactScalarState spExact;

      // Add second scalar random initial state generator
      spExact = spGen->addScalarEquation<Equations::ExactScalarState>();
      spExact->setIdentity(PhysicalNames::TEMPERATURE);
      spExact->setStateType(Equations::ExactScalarState::SINE);
      spExact->setSineOptions(10, 5, 5, 2);

      // Add second scalar random initial state generator
      spExact = spGen->addScalarEquation<Equations::ExactScalarState>();
      spExact->setIdentity(PhysicalNames::STREAMFUNCTION);
      spExact->setStateType(Equations::ExactScalarState::SINE);
      spExact->setSineOptions(10, 5, 5, 2);

      // Add second scalar random initial state generator
      spExact = spGen->addScalarEquation<Equations::ExactScalarState>();
      spExact->setIdentity(PhysicalNames::VORTICITY);
      spExact->setStateType(Equations::ExactScalarState::SINE);
      spExact->setSineOptions(10, 5, 5, 2);

      // Add second scalar random initial state generator
      spExact = spGen->addScalarEquation<Equations::ExactScalarState>();
      spExact->setIdentity(PhysicalNames::PRESSURE);
      spExact->setStateType(Equations::ExactScalarState::SINE);
      spExact->setSineOptions(10, 5, 5, 2);

      // Add output file
      IoVariable::SharedStateFileWriter spOut(new IoVariable::StateFileWriter(SchemeType::type(), SchemeType::isRegular()));
      spOut->expect(PhysicalNames::TEMPERATURE);
      spOut->expect(PhysicalNames::STREAMFUNCTION);
      spOut->expect(PhysicalNames::VORTICITY);
      spOut->expect(PhysicalNames::PRESSURE);
      spGen->addOutputFile(spOut);
   }

   void TestTTTModel::addVisualizers(SharedVisualizationGenerator spVis)
   {
      // Shared pointer to basic field visualizer
      Equations::SharedFieldVisualizer spField;

      // Add first field visualization
      spField = spVis->addScalarEquation<Equations::FieldVisualizer>();
      spField->setFields(true, false);
      spField->setIdentity(PhysicalNames::TEMPERATURE);

      // Add second field visualization
      spField = spVis->addScalarEquation<Equations::FieldVisualizer>();
      spField->setFields(true, false);
      spField->setIdentity(PhysicalNames::STREAMFUNCTION);

      // Add second field visualization
      spField = spVis->addScalarEquation<Equations::FieldVisualizer>();
      spField->setFields(true, false);
      spField->setIdentity(PhysicalNames::VORTICITY);

      // Add output file
      IoVariable::SharedVisualizationFileWriter spOut(new IoVariable::VisualizationFileWriter(SchemeType::type()));
      spOut->expect(PhysicalNames::TEMPERATURE);
      spOut->expect(PhysicalNames::STREAMFUNCTION);
      spOut->expect(PhysicalNames::VORTICITY);
      spVis->addOutputFile(spOut);
   }

   void TestTTTModel::setVisualizationState(SharedVisualizationGenerator spVis)
   {
      // Create and add initial state file to IO
      IoVariable::SharedStateFileReader spIn(new IoVariable::StateFileReader("4Visu", SchemeType::type(), SchemeType::isRegular()));

      // Set expected fields
      spIn->expect(PhysicalNames::TEMPERATURE);
      spIn->expect(PhysicalNames::STREAMFUNCTION);
      spIn->expect(PhysicalNames::VORTICITY);

      // Set simulation state
      spVis->setInitialState(spIn);
   }

   void TestTTTModel::addAsciiOutputFiles(SharedSimulation spSim)
   {
      // Add ASCII output file
      //pSim->addOutputFile(AN_ASCIIFILE);
      
      // Add ASCII output file
      //pSim->addOutputFile(AN_ASCIIFILE);
   }

   void TestTTTModel::addHdf5OutputFiles(SharedSimulation spSim)
   {
      // Field IDs iterator
      std::vector<GeoMHDiSCC::PhysicalNames::Id>::const_iterator  it;
      std::vector<GeoMHDiSCC::PhysicalNames::Id> ids = TestTTTModel::fieldIds();

      // Create and add state file to IO
      IoVariable::SharedStateFileWriter spState(new IoVariable::StateFileWriter(SchemeType::type(), SchemeType::isRegular()));
      for(it = ids.begin(); it != ids.end(); ++it)
      {
         spState->expect(*it);
      }
      spSim->addOutputFile(spState);
   }

   SharedSimulationBoundary TestTTTModel::createBoundary(const std::map<std::string,int>& bcIds)
   {
      // Create shared simulation boundary
      SharedSimulationBoundary  spBcs(new SimulationBoundary());

      // Storage for the dimension ID
      Dimensions::Simulation::Id dimId;

      // Create equation and field keys
      SpectralFieldId eqId;
      SpectralFieldId fieldId;

      // First scalar equation Dirichlet boundary conditions
      PhysicalNames::Id fieldName = PhysicalNames::TEMPERATURE;
      eqId = std::make_pair(fieldName, FieldComponents::Spectral::SCALAR);
      spBcs->initStorage(eqId);
      dimId = Dimensions::Simulation::SIM1D;
      fieldId = std::make_pair(fieldName, FieldComponents::Spectral::SCALAR);
      spBcs->initBcStorage(eqId, fieldId, dimId);
      spBcs->addBc(eqId, fieldId, dimId, Spectral::BoundaryConditions::VALUE, Spectral::IBoundary::LEFT);
      spBcs->addBc(eqId, fieldId, dimId, Spectral::BoundaryConditions::VALUE, Spectral::IBoundary::RIGHT);
      //spBcs->addBc(eqId, fieldId, dimId, Spectral::BoundaryConditions::SECOND_DERIVATIVE, Spectral::IBoundary::LEFT);
      //spBcs->addBc(eqId, fieldId, dimId, Spectral::BoundaryConditions::SECOND_DERIVATIVE, Spectral::IBoundary::RIGHT);
      dimId = Dimensions::Simulation::SIM3D;
      spBcs->initBcStorage(eqId, fieldId, dimId);
      spBcs->addBc(eqId, fieldId, dimId, Spectral::BoundaryConditions::VALUE, Spectral::IBoundary::LEFT);
      spBcs->addBc(eqId, fieldId, dimId, Spectral::BoundaryConditions::VALUE, Spectral::IBoundary::RIGHT);
      //spBcs->addBc(eqId, fieldId, dimId, Spectral::BoundaryConditions::SECOND_DERIVATIVE, Spectral::IBoundary::LEFT);
      //spBcs->addBc(eqId, fieldId, dimId, Spectral::BoundaryConditions::SECOND_DERIVATIVE, Spectral::IBoundary::RIGHT);

      // Second scalar equation Dirichlet boundary conditions
      fieldName = PhysicalNames::STREAMFUNCTION;
      eqId = std::make_pair(fieldName, FieldComponents::Spectral::SCALAR);
      spBcs->initStorage(eqId);
      dimId = Dimensions::Simulation::SIM1D;
      fieldId = std::make_pair(fieldName, FieldComponents::Spectral::SCALAR);
      spBcs->initBcStorage(eqId, fieldId, dimId);
      spBcs->addBc(eqId, fieldId, dimId, Spectral::BoundaryConditions::VALUE, Spectral::IBoundary::LEFT);
      spBcs->addBc(eqId, fieldId, dimId, Spectral::BoundaryConditions::VALUE, Spectral::IBoundary::RIGHT);
      spBcs->addBc(eqId, fieldId, dimId, Spectral::BoundaryConditions::SECOND_DERIVATIVE, Spectral::IBoundary::LEFT);
      spBcs->addBc(eqId, fieldId, dimId, Spectral::BoundaryConditions::SECOND_DERIVATIVE, Spectral::IBoundary::RIGHT);
      dimId = Dimensions::Simulation::SIM3D;
      spBcs->initBcStorage(eqId, fieldId, dimId);
      spBcs->addBc(eqId, fieldId, dimId, Spectral::BoundaryConditions::VALUE, Spectral::IBoundary::LEFT);
      spBcs->addBc(eqId, fieldId, dimId, Spectral::BoundaryConditions::VALUE, Spectral::IBoundary::RIGHT);
      //spBcs->addBc(eqId, fieldId, dimId, Spectral::BoundaryConditions::SECOND_DERIVATIVE, Spectral::IBoundary::LEFT);
      //spBcs->addBc(eqId, fieldId, dimId, Spectral::BoundaryConditions::SECOND_DERIVATIVE, Spectral::IBoundary::RIGHT);

      // Second scalar equation Dirichlet boundary conditions
      fieldName = PhysicalNames::VORTICITY;
      eqId = std::make_pair(fieldName, FieldComponents::Spectral::SCALAR);
      spBcs->initStorage(eqId);
      dimId = Dimensions::Simulation::SIM1D;
      fieldId = std::make_pair(fieldName, FieldComponents::Spectral::SCALAR);
      spBcs->initBcStorage(eqId, fieldId, dimId);
      spBcs->addBc(eqId, fieldId, dimId, Spectral::BoundaryConditions::VALUE, Spectral::IBoundary::LEFT);
      spBcs->addBc(eqId, fieldId, dimId, Spectral::BoundaryConditions::VALUE, Spectral::IBoundary::RIGHT);
      fieldName = PhysicalNames::PRESSURE;
      eqId = std::make_pair(fieldName, FieldComponents::Spectral::SCALAR);
      spBcs->initStorage(eqId);
      dimId = Dimensions::Simulation::SIM1D;
      fieldId = std::make_pair(PhysicalNames::VORTICITY, FieldComponents::Spectral::SCALAR);
      spBcs->initBcStorage(eqId, fieldId, dimId);
      spBcs->addBc(eqId, fieldId, dimId, Spectral::BoundaryConditions::SECOND_DERIVATIVE, Spectral::IBoundary::LEFT);
      spBcs->addBc(eqId, fieldId, dimId, Spectral::BoundaryConditions::SECOND_DERIVATIVE, Spectral::IBoundary::RIGHT);

      return spBcs;
   }

   void TestTTTModel::setInitialState(SharedSimulation spSim)
   {
      // Field IDs iterator
      std::vector<GeoMHDiSCC::PhysicalNames::Id>::const_iterator  it;
      std::vector<GeoMHDiSCC::PhysicalNames::Id> ids = TestTTTModel::fieldIds();

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
