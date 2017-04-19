/** 
 * @file PhysicalModel.cpp
 * @brief Source of the test model for the FFF scheme
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

/// Define small macros allowing to convert to string
#define MAKE_STR_X( _P ) # _P
#define MAKE_STR( _P ) MAKE_STR_X( _P )

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include MAKE_STR( QUICC_MODEL_PATH/Test/FFF/Linear/PhysicalModel.hpp )

// Project includes
//
#include MAKE_STR( QUICC_MODEL_PATH/Test/LinearScalar.hpp )
#include MAKE_STR( QUICC_MODEL_PATH/Test/NonlinearScalar.hpp )
#include "Enums/FieldIds.hpp"
#include "IoVariable/StateFileReader.hpp"
#include "IoVariable/StateFileWriter.hpp"
#include "IoVariable/VisualizationFileWriter.hpp"
#include "IoTools/IdToHuman.hpp"
#include "Generator/States/RandomScalarState.hpp"
#include "Generator/States/CartesianExactScalarState.hpp"
#include "Generator/Visualizers/ScalarFieldVisualizer.hpp"
#include "Model/PhysicalModelBase.hpp"

namespace QuICC {

namespace Model {

namespace Test {

namespace FFF {

namespace Linear {

   const std::string PhysicalModel::PYMODULE = "test_fff";

   const std::string PhysicalModel::PYCLASS = "TestFFF";

   void PhysicalModel::addEquations(SharedSimulation spSim)
   {
      Equations::Test::SharedLinearScalar   spLin;

      // Add first scalar test equation
      spLin = spSim->addScalarEquation<Equations::Test::LinearScalar>();
      spLin->setIdentity(PhysicalNames::STREAMFUNCTION);

      // Add second scalar test equation
      spLin = spSim->addScalarEquation<Equations::Test::LinearScalar>();
      spLin->setIdentity(PhysicalNames::VELOCITYZ);

      // Add third scalar test equation
      spLin = spSim->addScalarEquation<Equations::Test::LinearScalar>();
      spLin->setIdentity(PhysicalNames::TEMPERATURE);
   }

   void PhysicalModel::addStates(SharedStateGenerator spGen)
   {
      // Shared pointer to equation
      Equations::SharedRandomScalarState spRand;
      Equations::SharedCartesianExactScalarState spExact;

      // Add first scalar initial state generator
      spExact = spGen->addScalarEquation<Equations::CartesianExactScalarState>();
      spExact->setIdentity(PhysicalNames::STREAMFUNCTION);
      spExact->setStateType(Equations::CartesianExactStateIds::SINSINCOS);
      spExact->setModeOptions(1e0, 1.0, 1e0, 1.0, 1e0, 1.0);

      // Add first scalar initial state generator
      spExact = spGen->addScalarEquation<Equations::CartesianExactScalarState>();
      spExact->setIdentity(PhysicalNames::VELOCITYZ);
      spExact->setStateType(Equations::CartesianExactStateIds::COSSINCOS);
      spExact->setModeOptions(1e0, 1.0, 1e0, 1.0, 1e0, 1.0);

      // Add first scalar initial state generator
      spExact = spGen->addScalarEquation<Equations::CartesianExactScalarState>();
      spExact->setIdentity(PhysicalNames::TEMPERATURE);
      spExact->setStateType(Equations::CartesianExactStateIds::SINCOSCOS);
      spExact->setModeOptions(1e0, 1.0, 1e0, 1.0, 1e0, 1.0);

      // Add output file
      IoVariable::SharedStateFileWriter spOut(new IoVariable::StateFileWriter(SchemeType::type(), SchemeType::isRegular()));
      spOut->expect(PhysicalNames::STREAMFUNCTION);
      spOut->expect(PhysicalNames::VELOCITYZ);
      spOut->expect(PhysicalNames::TEMPERATURE);
      spGen->addHdf5OutputFile(spOut);
   }

   void PhysicalModel::addVisualizers(SharedVisualizationGenerator spVis)
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

   void PhysicalModel::setVisualizationState(SharedVisualizationGenerator spVis)
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

   void PhysicalModel::addAsciiOutputFiles(SharedSimulation spSim)
   {
      // Add ASCII output file
      //pSim->addOutputFile(AN_ASCIIFILE);
      
      // Add ASCII output file
      //pSim->addOutputFile(AN_ASCIIFILE);
   }

   void PhysicalModel::addHdf5OutputFiles(SharedSimulation spSim)
   {
      // Field IDs iterator
      std::vector<PhysicalNames::Id>::const_iterator  it;
      std::vector<PhysicalNames::Id> ids = PhysicalModelBase::fieldIds();

      // Create and add state file to IO
      IoVariable::SharedStateFileWriter spState(new IoVariable::StateFileWriter(SchemeType::type(), SchemeType::isRegular()));
      for(it = ids.begin(); it != ids.end(); ++it)
      {
         spState->expect(*it);
      }
      spSim->addHdf5OutputFile(spState);
   }

   void PhysicalModel::setInitialState(SharedSimulation spSim)
   {
      // Field IDs iterator
      std::vector<PhysicalNames::Id>::const_iterator  it;
      std::vector<PhysicalNames::Id> ids = PhysicalModelBase::fieldIds();

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
}
}
}
}
