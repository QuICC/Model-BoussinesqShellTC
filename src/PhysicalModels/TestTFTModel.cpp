/** \file TestTFTModel.cpp
 *  \brief Source of the test model for the TFT scheme
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "PhysicalModels/TestTFTModel.hpp"

// Project includes
//
#include "Enums/FieldIds.hpp"
#include "IoVariable/StateFileReader.hpp"
#include "IoVariable/StateFileWriter.hpp"
#include "IoVariable/VisualizationFileWriter.hpp"
#include "IoTools/IdToHuman.hpp"
#include "Equations/Tests/TestTFTScalarOne.hpp"
#include "Equations/Tests/TestTFTScalarTwo.hpp"
#include "Equations/Tests/TestTFTScalarThree.hpp"

namespace GeoMHDiSCC {

   std::vector<PhysicalNames::Id> TestTFTModel::fieldIds()
   {
      // Create storage
      std::vector<PhysicalNames::Id> ids;

      // Add first scalar
      ids.push_back(PhysicalNames::TEMPERATURE);

      // Add second scalar
      ids.push_back(PhysicalNames::STREAMFUNCTION);

      // Add third scalar
      ids.push_back(PhysicalNames::VELOCITYZ);

      return ids;
   }

   std::vector<NonDimensional::Id> TestTFTModel::paramIds()
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

   std::vector<bool> TestTFTModel::isPeriodicBox()
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

   void TestTFTModel::addEquations(SharedSimulation spSim)
   {
      // Add first scalar test equation
      spSim->addScalarEquation<Equations::TestTFTScalarOne>();

      // Add second scalar test equation
      spSim->addScalarEquation<Equations::TestTFTScalarTwo>();

      // Add third scalar test equation
      spSim->addScalarEquation<Equations::TestTFTScalarThree>();
   }

   void TestTFTModel::addAsciiOutputFiles(SharedSimulation spSim)
   {
      // Add ASCII output file
      //pSim->addOutputFile(AN_ASCIIFILE);
      
      // Add ASCII output file
      //pSim->addOutputFile(AN_ASCIIFILE);
   }

   void TestTFTModel::addHdf5OutputFiles(SharedSimulation spSim)
   {
      // Field IDs iterator
      std::vector<GeoMHDiSCC::PhysicalNames::Id>::const_iterator  it;
      std::vector<GeoMHDiSCC::PhysicalNames::Id> ids = TestTFTModel::fieldIds();

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

   SharedSimulationBoundary TestTFTModel::createBoundary(const std::map<std::string,int>& bcIds)
   {
      // Create shared simulation boundary
      SharedSimulationBoundary  spBcs(new SimulationBoundary());

      // Storage for the dimension ID
      Dimensions::Simulation::Id dimId;

      // Create equation and field keys
      SpectralFieldId eqId;
      SpectralFieldId fieldId;

      // First scalar equation Dirichlet boundary conditions
      eqId = std::make_pair(PhysicalNames::TEMPERATURE, FieldComponents::Spectral::SCALAR);
      spBcs->initStorage(eqId);
      dimId = Dimensions::Simulation::SIM1D;
      fieldId = std::make_pair(PhysicalNames::TEMPERATURE, FieldComponents::Spectral::SCALAR);
      spBcs->initBcStorage(eqId, fieldId, dimId);
      spBcs->addBc(eqId, fieldId, dimId, Spectral::BoundaryConditions::VALUE, Spectral::IBoundary::LEFT);
      spBcs->addBc(eqId, fieldId, dimId, Spectral::BoundaryConditions::VALUE, Spectral::IBoundary::RIGHT);
      dimId = Dimensions::Simulation::SIM3D;
      spBcs->initBcStorage(eqId, fieldId, dimId);
      spBcs->addBc(eqId, fieldId, dimId, Spectral::BoundaryConditions::VALUE, Spectral::IBoundary::LEFT);
      spBcs->addBc(eqId, fieldId, dimId, Spectral::BoundaryConditions::VALUE, Spectral::IBoundary::RIGHT);

      // Second scalar equation Dirichlet boundary conditions
      eqId = std::make_pair(PhysicalNames::STREAMFUNCTION, FieldComponents::Spectral::SCALAR);
      spBcs->initStorage(eqId);
      dimId = Dimensions::Simulation::SIM1D;
      fieldId = std::make_pair(PhysicalNames::STREAMFUNCTION, FieldComponents::Spectral::SCALAR);
      spBcs->initBcStorage(eqId, fieldId, dimId);
      spBcs->addBc(eqId, fieldId, dimId, Spectral::BoundaryConditions::VALUE, Spectral::IBoundary::LEFT);
      spBcs->addBc(eqId, fieldId, dimId, Spectral::BoundaryConditions::VALUE, Spectral::IBoundary::RIGHT);
      dimId = Dimensions::Simulation::SIM3D;
      spBcs->initBcStorage(eqId, fieldId, dimId);
      spBcs->addBc(eqId, fieldId, dimId, Spectral::BoundaryConditions::VALUE, Spectral::IBoundary::LEFT);
      spBcs->addBc(eqId, fieldId, dimId, Spectral::BoundaryConditions::VALUE, Spectral::IBoundary::RIGHT);

      // Third scalar equation Dirichlet boundary conditions
      eqId = std::make_pair(PhysicalNames::VELOCITYZ, FieldComponents::Spectral::SCALAR);
      spBcs->initStorage(eqId);
      dimId = Dimensions::Simulation::SIM1D;
      fieldId = std::make_pair(PhysicalNames::VELOCITYZ, FieldComponents::Spectral::SCALAR);
      spBcs->initBcStorage(eqId, fieldId, dimId);
      spBcs->addBc(eqId, fieldId, dimId, Spectral::BoundaryConditions::VALUE, Spectral::IBoundary::LEFT);
      spBcs->addBc(eqId, fieldId, dimId, Spectral::BoundaryConditions::VALUE, Spectral::IBoundary::RIGHT);
      dimId = Dimensions::Simulation::SIM3D;
      spBcs->initBcStorage(eqId, fieldId, dimId);
      spBcs->addBc(eqId, fieldId, dimId, Spectral::BoundaryConditions::VALUE, Spectral::IBoundary::LEFT);
      spBcs->addBc(eqId, fieldId, dimId, Spectral::BoundaryConditions::VALUE, Spectral::IBoundary::RIGHT);

      return spBcs;
   }

   void TestTFTModel::setInitialState(SharedSimulation spSim)
   {
      // Field IDs iterator
      std::vector<GeoMHDiSCC::PhysicalNames::Id>::const_iterator  it;
      std::vector<GeoMHDiSCC::PhysicalNames::Id> ids = TestTFTModel::fieldIds();

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
