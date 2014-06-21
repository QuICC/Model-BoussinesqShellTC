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
#include "BoundaryCondition/BoundaryCondition.hpp"
#include "IoVariable/StateFileReader.hpp"
#include "IoVariable/StateFileWriter.hpp"
#include "IoVariable/VisualizationFileWriter.hpp"
#include "IoTools/IdToHuman.hpp"
#include "Equations/Tests/TestLinearScalar.hpp"
#include "Generator/States/RandomScalarState.hpp"
#include "Generator/States/CartesianExactScalarState.hpp"
#include "Generator/Visualizers/ScalarFieldVisualizer.hpp"
#include "PhysicalModels/PhysicalModelBase.hpp"

namespace GeoMHDiSCC {

   const std::string TestTTTModel::PYMODULE = "test_ttt_model";

   const std::string TestTTTModel::PYCLASS = "TestTTTModel";

   void TestTTTModel::addEquations(SharedSimulation spSim)
   {
      Equations::SharedTestLinearScalar   spLin;

      // Add first scalar test equation
      spLin = spSim->addScalarEquation<Equations::TestLinearScalar>();
      spLin->setIdentity(PhysicalNames::STREAMFUNCTION);

      // Add second scalar test equation
      spLin = spSim->addScalarEquation<Equations::TestLinearScalar>();
      spLin->setIdentity(PhysicalNames::VELOCITYZ);

      // Add third scalar test equation
      spLin = spSim->addScalarEquation<Equations::TestLinearScalar>();
      spLin->setIdentity(PhysicalNames::TEMPERATURE);
   }

   void TestTTTModel::addStates(SharedStateGenerator spGen)
   {
      // Shared pointer to equation
      Equations::SharedRandomScalarState spRand;
      Equations::SharedCartesianExactScalarState spExact;

      // Add first scalar initial state generator
      spExact = spGen->addScalarEquation<Equations::CartesianExactScalarState>();
      spExact->setIdentity(PhysicalNames::STREAMFUNCTION);
      spExact->setStateType(Equations::CartesianExactScalarState::POLYPOLYPOLY);
      spExact->setModeOptions(1e0, 1.0, 1e0, 1.0, 1e0, 1.0);

      // Add first scalar initial state generator
      spExact = spGen->addScalarEquation<Equations::CartesianExactScalarState>();
      spExact->setIdentity(PhysicalNames::VELOCITYZ);
      spExact->setStateType(Equations::CartesianExactScalarState::POLYPOLYPOLY);
      spExact->setModeOptions(1e0, 1.0, 1e0, 1.0, 1e0, 1.0);

      // Add first scalar initial state generator
      spExact = spGen->addScalarEquation<Equations::CartesianExactScalarState>();
      spExact->setIdentity(PhysicalNames::TEMPERATURE);
      spExact->setStateType(Equations::CartesianExactScalarState::POLYPOLYPOLY);
      spExact->setModeOptions(1e0, 1.0, 1e0, 1.0, 1e0, 1.0);

      // Add output file
      IoVariable::SharedStateFileWriter spOut(new IoVariable::StateFileWriter(SchemeType::type(), SchemeType::isRegular()));
      spOut->expect(PhysicalNames::STREAMFUNCTION);
      spOut->expect(PhysicalNames::VELOCITYZ);
      spOut->expect(PhysicalNames::TEMPERATURE);
      spGen->addHdf5OutputFile(spOut);
   }

   void TestTTTModel::addVisualizers(SharedVisualizationGenerator spVis)
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
      spVis->addHdf5OutputFile(spOut);
   }

   void TestTTTModel::setVisualizationState(SharedVisualizationGenerator spVis)
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
      std::vector<GeoMHDiSCC::PhysicalNames::Id> ids = PhysicalModelBase::fieldIds();

      // Create and add state file to IO
      IoVariable::SharedStateFileWriter spState(new IoVariable::StateFileWriter(SchemeType::type(), SchemeType::isRegular()));
      for(it = ids.begin(); it != ids.end(); ++it)
      {
         spState->expect(*it);
      }
      spSim->addHdf5OutputFile(spState);
   }

   void TestTTTModel::setInitialState(SharedSimulation spSim)
   {
      // Field IDs iterator
      std::vector<GeoMHDiSCC::PhysicalNames::Id>::const_iterator  it;
      std::vector<GeoMHDiSCC::PhysicalNames::Id> ids = PhysicalModelBase::fieldIds();

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
