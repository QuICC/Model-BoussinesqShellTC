/** 
 * @file TestFFFModel.cpp
 * @brief Source of the test model for the FFF scheme
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
#include "PhysicalModels/TestFFFModel.hpp"

// Project includes
//
#include "Enums/FieldIds.hpp"
#include "BoundaryCondition/BoundaryCondition.hpp"
#include "IoVariable/StateFileReader.hpp"
#include "IoVariable/StateFileWriter.hpp"
#include "IoVariable/VisualizationFileWriter.hpp"
#include "IoTools/IdToHuman.hpp"
#include "Equations/Tests/TestFFFDiffusion2D.hpp"
#include "Equations/Tests/TestFFFDiffusion3D.hpp"
#include "Generator/States/RandomScalarState.hpp"
#include "Generator/States/ExactScalarState.hpp"
#include "Generator/Visualizers/ScalarFieldVisualizer.hpp"

namespace GeoMHDiSCC {

   const std::string TestFFFModel::PYNAME = "test_fff_model.py";

   std::vector<PhysicalNames::Id> TestFFFModel::fieldIds()
   {
      // Create storage
      std::vector<PhysicalNames::Id> ids;

      // Add second scalar
      ids.push_back(PhysicalNames::STREAMFUNCTION);

      // Add third scalar
      ids.push_back(PhysicalNames::VELOCITYZ);

      // Add first scalar
      ids.push_back(PhysicalNames::TEMPERATURE);

      return ids;
   }

   std::vector<NonDimensional::Id> TestFFFModel::paramIds()
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

   std::vector<bool> TestFFFModel::isPeriodicBox()
   {
      std::vector<bool> box;

      // X direction is not periodic box
      box.push_back(true);

      // Y direction is periodic box
      box.push_back(true);

      // Z direction is not periodic box
      box.push_back(true);

      return box;
   }

   void TestFFFModel::addEquations(SharedSimulation spSim)
   {
      Equations::SharedTestFFFDiffusion2D   spD2D;

      // Add first scalar test equation
      spD2D = spSim->addScalarEquation<Equations::TestFFFDiffusion2D>(TestFFFModel::PYNAME);
      spD2D->setIdentity(PhysicalNames::STREAMFUNCTION);

      Equations::SharedTestFFFDiffusion3D   spD3D;

      // Add second scalar test equation
      spD3D = spSim->addScalarEquation<Equations::TestFFFDiffusion3D>(TestFFFModel::PYNAME);
      spD3D->setIdentity(PhysicalNames::VELOCITYZ);

      // Add third scalar test equation
      spD2D = spSim->addScalarEquation<Equations::TestFFFDiffusion2D>(TestFFFModel::PYNAME);
      spD2D->setIdentity(PhysicalNames::TEMPERATURE);
   }

   void TestFFFModel::addStates(SharedStateGenerator spGen)
   {
      // Shared pointer to equation
      Equations::SharedRandomScalarState spRand;
      Equations::SharedExactScalarState spExact;

      // Add first scalar initial state generator
      spExact = spGen->addScalarEquation<Equations::ExactScalarState>();
      spExact->setIdentity(PhysicalNames::STREAMFUNCTION);
      spExact->setStateType(Equations::ExactScalarState::SINECOSINE);
      spExact->setSineOptions(10, 5, 5, 2);

      // Add first scalar initial state generator
      spExact = spGen->addScalarEquation<Equations::ExactScalarState>();
      spExact->setIdentity(PhysicalNames::VELOCITYZ);
      spExact->setStateType(Equations::ExactScalarState::SINECOSINE);
      spExact->setSineOptions(10, 5, 5, 2);

      // Add first scalar initial state generator
      spExact = spGen->addScalarEquation<Equations::ExactScalarState>();
      spExact->setIdentity(PhysicalNames::TEMPERATURE);
      spExact->setStateType(Equations::ExactScalarState::SINECOSINE);
      spExact->setSineOptions(10, 5, 5, 2);

      // Add output file
      IoVariable::SharedStateFileWriter spOut(new IoVariable::StateFileWriter(SchemeType::type(), SchemeType::isRegular()));
      spOut->expect(PhysicalNames::STREAMFUNCTION);
      spOut->expect(PhysicalNames::VELOCITYZ);
      spOut->expect(PhysicalNames::TEMPERATURE);
      spGen->addOutputFile(spOut);
   }

   void TestFFFModel::addVisualizers(SharedVisualizationGenerator spVis)
   {
      // Shared pointer to basic field visualizer
      Equations::SharedScalarFieldVisualizer spField;

      // Add first field visualization
      spField = spVis->addScalarEquation<Equations::ScalarFieldVisualizer>();
      spField->setFields(true, false);
      spField->setIdentity(PhysicalNames::STREAMFUNCTION);

      // Add second field visualization
      spField = spVis->addScalarEquation<Equations::ScalarFieldVisualizer>();
      spField->setFields(true, false);
      spField->setIdentity(PhysicalNames::VELOCITYZ);

      // Add third field visualization
      spField = spVis->addScalarEquation<Equations::ScalarFieldVisualizer>();
      spField->setFields(true, false);
      spField->setIdentity(PhysicalNames::TEMPERATURE);

      // Add output file
      IoVariable::SharedVisualizationFileWriter spOut(new IoVariable::VisualizationFileWriter(SchemeType::type()));
      spOut->expect(PhysicalNames::STREAMFUNCTION);
      spOut->expect(PhysicalNames::VELOCITYZ);
      spOut->expect(PhysicalNames::TEMPERATURE);
      spVis->addOutputFile(spOut);
   }

   void TestFFFModel::setVisualizationState(SharedVisualizationGenerator spVis)
   {
      // Create input state file for visualisation
      IoVariable::SharedStateFileReader spIn(new IoVariable::StateFileReader("4Visu", SchemeType::type(), SchemeType::isRegular()));

      // Set expected fields
      spIn->expect(PhysicalNames::STREAMFUNCTION);
      spIn->expect(PhysicalNames::TEMPERATURE);
      spIn->expect(PhysicalNames::VELOCITYZ);

      // Set state for visualization generator
      spVis->setInitialState(spIn);
   }

   void TestFFFModel::addAsciiOutputFiles(SharedSimulation spSim)
   {
      // Add ASCII output file
      //pSim->addOutputFile(AN_ASCIIFILE);
      
      // Add ASCII output file
      //pSim->addOutputFile(AN_ASCIIFILE);
   }

   void TestFFFModel::addHdf5OutputFiles(SharedSimulation spSim)
   {
      // Field IDs iterator
      std::vector<GeoMHDiSCC::PhysicalNames::Id>::const_iterator  it;
      std::vector<GeoMHDiSCC::PhysicalNames::Id> ids = TestFFFModel::fieldIds();

      // Create and add state file to IO
      IoVariable::SharedStateFileWriter spState(new IoVariable::StateFileWriter(SchemeType::type(), SchemeType::isRegular()));
      for(it = ids.begin(); it != ids.end(); ++it)
      {
         spState->expect(*it);
      }
      spSim->addOutputFile(spState);
   }

   SharedSimulationBoundary TestFFFModel::createBoundary(const std::map<std::string,int>& bcIds)
   {
      // Create shared simulation boundary
      SharedSimulationBoundary  spBcs(new SimulationBoundary());

      // Storage for the dimension ID
      Dimensions::Simulation::Id dimId;

      // Create equation and field keys
      SpectralFieldId eqId;
      SpectralFieldId fieldId;

      // Second scalar equation Dirichlet boundary conditions
      PhysicalNames::Id fieldName = PhysicalNames::STREAMFUNCTION;
      eqId = std::make_pair(fieldName, FieldComponents::Spectral::SCALAR);
      spBcs->initStorage(eqId);
      dimId = Dimensions::Simulation::SIM1D;
      fieldId = std::make_pair(fieldName, FieldComponents::Spectral::SCALAR);
      spBcs->initBcStorage(eqId, fieldId, dimId);
      spBcs->addBc(eqId, fieldId, dimId, Boundary::VALUE, Boundary::LEFT);
      spBcs->addBc(eqId, fieldId, dimId, Boundary::VALUE, Boundary::RIGHT);
      dimId = Dimensions::Simulation::SIM3D;
      spBcs->initBcStorage(eqId, fieldId, dimId);
      spBcs->addBc(eqId, fieldId, dimId, Boundary::VALUE, Boundary::LEFT);
      spBcs->addBc(eqId, fieldId, dimId, Boundary::VALUE, Boundary::RIGHT);

      // First scalar equation Dirichlet boundary conditions
      fieldName = PhysicalNames::TEMPERATURE;
      eqId = std::make_pair(fieldName, FieldComponents::Spectral::SCALAR);
      spBcs->initStorage(eqId);
      dimId = Dimensions::Simulation::SIM1D;
      fieldId = std::make_pair(fieldName, FieldComponents::Spectral::SCALAR);
      spBcs->initBcStorage(eqId, fieldId, dimId);
      spBcs->addBc(eqId, fieldId, dimId, Boundary::VALUE, Boundary::LEFT);
      spBcs->addBc(eqId, fieldId, dimId, Boundary::VALUE, Boundary::RIGHT);
      dimId = Dimensions::Simulation::SIM3D;
      spBcs->initBcStorage(eqId, fieldId, dimId);
      spBcs->addBc(eqId, fieldId, dimId, Boundary::VALUE, Boundary::LEFT);
      spBcs->addBc(eqId, fieldId, dimId, Boundary::VALUE, Boundary::RIGHT);

      // Second scalar equation Dirichlet boundary conditions
      fieldName = PhysicalNames::VELOCITYZ;
      eqId = std::make_pair(fieldName, FieldComponents::Spectral::SCALAR);
      spBcs->initStorage(eqId);
      dimId = Dimensions::Simulation::SIM1D;
      fieldId = std::make_pair(fieldName, FieldComponents::Spectral::SCALAR);
      spBcs->initBcStorage(eqId, fieldId, dimId);
      spBcs->addBc(eqId, fieldId, dimId, Boundary::VALUE, Boundary::LEFT);
      spBcs->addBc(eqId, fieldId, dimId, Boundary::VALUE, Boundary::RIGHT);
      dimId = Dimensions::Simulation::SIM3D;
      spBcs->initBcStorage(eqId, fieldId, dimId);
      spBcs->addBc(eqId, fieldId, dimId, Boundary::VALUE, Boundary::LEFT);
      spBcs->addBc(eqId, fieldId, dimId, Boundary::VALUE, Boundary::RIGHT);

      return spBcs;
   }

   void TestFFFModel::setInitialState(SharedSimulation spSim)
   {
      // Field IDs iterator
      std::vector<GeoMHDiSCC::PhysicalNames::Id>::const_iterator  it;
      std::vector<GeoMHDiSCC::PhysicalNames::Id> ids = TestFFFModel::fieldIds();

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
