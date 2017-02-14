/** 
 * @file PhysicalModel.cpp
 * @brief Source of the test model for the TFF scheme
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
#include MAKE_STR( QUICC_MODEL_PATH/Test/TFF/Linear/PhysicalModel.hpp )

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

namespace TFF {

namespace Linear {

   const std::string PhysicalModel::PYMODULE = "test_tffscheme";

   const std::string PhysicalModel::PYCLASS = "TestTFFScheme";

   void PhysicalModel::addEquations(SharedSimulation spSim)
   {
      // Create linear test problem
      if(true)
      {
         Equations::Test::SharedLinearScalar   spLin;

//         // Add scalar test equation
//         spLin = spSim->addScalarEquation<Equations::Test::LinearScalar>();
//         spLin->setIdentity(PhysicalNames::VELOCITYX);
//
//         // Add scalar test equation
//         spLin = spSim->addScalarEquation<Equations::Test::LinearScalar>();
//         spLin->setIdentity(PhysicalNames::VELOCITYY);
//
//         // Add scalar test equation
//         spLin = spSim->addScalarEquation<Equations::Test::LinearScalar>();
//         spLin->setIdentity(PhysicalNames::VELOCITYZ);

         // Add scalar test equation
         spLin = spSim->addScalarEquation<Equations::Test::LinearScalar>();
         spLin->setIdentity(PhysicalNames::TEMPERATURE);
      } else
      {
         Equations::Test::SharedNonlinearScalar   spNL;

//         // Add scalar test equation
//         spNL = spSim->addScalarEquation<Equations::Test::NonlinearScalar>();
//         spNL->setIdentity(PhysicalNames::VELOCITYX);
//
//         // Add scalar test equation
//         spNL = spSim->addScalarEquation<Equations::Test::NonlinearScalar>();
//         spNL->setIdentity(PhysicalNames::VELOCITYY);
//
//         // Add scalar test equation
//         spNL = spSim->addScalarEquation<Equations::Test::NonlinearScalar>();
//         spNL->setIdentity(PhysicalNames::VELOCITYZ);

         // Add scalar test equation
         spNL = spSim->addScalarEquation<Equations::Test::NonlinearScalar>();
         spNL->setIdentity(PhysicalNames::TEMPERATURE);
      }
   }

   void PhysicalModel::addStates(SharedStateGenerator spGen)
   {
      // Generate "exact" solutions (trigonometric or monomial)
      if(true)
      {
         // Shared pointer to equation
         Equations::SharedCartesianExactScalarState spExact;

         // Add scalar exact initial state generator
         spExact = spGen->addScalarEquation<Equations::CartesianExactScalarState>();
         spExact->setIdentity(PhysicalNames::VELOCITYX);
         spExact->setStateType(Equations::CartesianExactStateIds::POLYSINCOS);
         spExact->setModeOptions(1e0, 0.0, 1e0, 1.0, 1e0, 1.0);

         // Add scalar exact initial state generator
         spExact = spGen->addScalarEquation<Equations::CartesianExactScalarState>();
         spExact->setIdentity(PhysicalNames::VELOCITYY);
         spExact->setStateType(Equations::CartesianExactStateIds::POLYSINCOS);
         spExact->setModeOptions(-1e0, 1.0, -2e0, 1.0, 1e1, 2.0);

         // Add scalar exact initial state generator
         spExact = spGen->addScalarEquation<Equations::CartesianExactScalarState>();
         spExact->setIdentity(PhysicalNames::VELOCITYZ);
         spExact->setStateType(Equations::CartesianExactStateIds::POLYSINCOS);
         spExact->setModeOptions(1e1, 2.0, 5e0, 2.0, 2e1, 1.0);

         // Add scalar exact initial state generator
         spExact = spGen->addScalarEquation<Equations::CartesianExactScalarState>();
         spExact->setIdentity(PhysicalNames::TEMPERATURE);
         spExact->setStateType(Equations::CartesianExactStateIds::POLYSINSIN);
         spExact->setModeOptions(-1e2, 10.0, 3e0, 10.0, -3e1, 10.0);

      // Generate random spectrum
      } else
      {
         // Shared pointer to random initial state equation
         Equations::SharedRandomScalarState spRand;

         // Add scalar random initial state generator 
         spRand = spGen->addScalarEquation<Equations::RandomScalarState>();
         spRand->setIdentity(PhysicalNames::VELOCITYX);
         spRand->setSpectrum(-0.001, 0.001, 1e4, 1e4, 1e4);

         // Add scalar random initial state generator 
         spRand = spGen->addScalarEquation<Equations::RandomScalarState>();
         spRand->setIdentity(PhysicalNames::VELOCITYY);
         spRand->setSpectrum(-0.001, 0.001, 1e4, 1e4, 1e4);

         // Add scalar random initial state generator 
         spRand = spGen->addScalarEquation<Equations::RandomScalarState>();
         spRand->setIdentity(PhysicalNames::VELOCITYZ);
         spRand->setSpectrum(-0.001, 0.001, 1e4, 1e4, 1e4);

         // Add scalar random initial state generator
         spRand = spGen->addScalarEquation<Equations::RandomScalarState>();
         spRand->setIdentity(PhysicalNames::TEMPERATURE);
         spRand->setSpectrum(-0.001, 0.001, 1e4, 1e4, 1e4);
      }

      // Add output file
      IoVariable::SharedStateFileWriter spOut(new IoVariable::StateFileWriter(SchemeType::type(), SchemeType::isRegular()));
      spOut->expect(PhysicalNames::VELOCITYX);
      spOut->expect(PhysicalNames::VELOCITYY);
      spOut->expect(PhysicalNames::VELOCITYZ);
      spOut->expect(PhysicalNames::TEMPERATURE);
      spGen->addHdf5OutputFile(spOut);
   }

   void PhysicalModel::addVisualizers(SharedVisualizationGenerator spVis)
   {
      // Shared pointer to basic field visualizer
      Equations::SharedScalarFieldVisualizer spField;

      // Add scalar field visualization
      spField = spVis->addScalarEquation<Equations::ScalarFieldVisualizer>();
      spField->setFields(true, true);
      spField->setIdentity(PhysicalNames::VELOCITYX);

      // Add scalar field visualization
      spField = spVis->addScalarEquation<Equations::ScalarFieldVisualizer>();
      spField->setFields(true, true);
      spField->setIdentity(PhysicalNames::VELOCITYY);

      // Add scalar field visualization
      spField = spVis->addScalarEquation<Equations::ScalarFieldVisualizer>();
      spField->setFields(true, true);
      spField->setIdentity(PhysicalNames::VELOCITYZ);

      // Add scalar field visualization
      spField = spVis->addScalarEquation<Equations::ScalarFieldVisualizer>();
      spField->setFields(true, true);
      spField->setIdentity(PhysicalNames::TEMPERATURE);

      // Add output file
      IoVariable::SharedVisualizationFileWriter spOut(new IoVariable::VisualizationFileWriter(SchemeType::type()));
      spOut->expect(PhysicalNames::VELOCITYX);
      spOut->expect(PhysicalNames::VELOCITYY);
      spOut->expect(PhysicalNames::VELOCITYZ);
      spOut->expect(PhysicalNames::TEMPERATURE);
      spVis->addHdf5OutputFile(spOut);
   }

   void PhysicalModel::setVisualizationState(SharedVisualizationGenerator spVis)
   {
      // Create input state file for visualisation
      IoVariable::SharedStateFileReader spIn(new IoVariable::StateFileReader("4Visu", SchemeType::type(), SchemeType::isRegular()));

      // Set expected fields
      spIn->expect(PhysicalNames::VELOCITYX);
      spIn->expect(PhysicalNames::VELOCITYY);
      spIn->expect(PhysicalNames::VELOCITYZ);
      spIn->expect(PhysicalNames::TEMPERATURE);

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

   void PhysicalModel::addStatsOutputFiles(SharedSimulation spSim)
   {
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
