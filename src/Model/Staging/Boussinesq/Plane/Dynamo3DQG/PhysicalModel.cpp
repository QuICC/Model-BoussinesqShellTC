/** 
 * @file PhysicalModel.cpp
 * @brief Source of the Boussinesq F-plane 3DQG physical model
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
#include MAKE_STR( QUICC_MODEL_PATH/Boussinesq/Plane/Dynamo3DQG/PhysicalModel.hpp )

// Project includes
//
#include MAKE_STR( QUICC_MODEL_PATH/Boussinesq/Plane/Dynamo3DQG/Streamfunction.hpp )
#include MAKE_STR( QUICC_MODEL_PATH/Boussinesq/Plane/Dynamo3DQG/VelocityZ.hpp )
#include MAKE_STR( QUICC_MODEL_PATH/Boussinesq/Plane/Dynamo3DQG/Transport.hpp )
#include MAKE_STR( QUICC_MODEL_PATH/Boussinesq/Plane/Dynamo3DQG/VorticityZ.hpp )
#include MAKE_STR( QUICC_MODEL_PATH/Boussinesq/Plane/Dynamo3DQG/Pressure.hpp )
#include MAKE_STR( QUICC_MODEL_PATH/Boussinesq/Plane/Dynamo3DQG/MeanHeat.hpp )
#include MAKE_STR( QUICC_MODEL_PATH/Boussinesq/Plane/Dynamo3DQG/Emfx.hpp )
#include MAKE_STR( QUICC_MODEL_PATH/Boussinesq/Plane/Dynamo3DQG/Emfy.hpp )
#include MAKE_STR( QUICC_MODEL_PATH/Boussinesq/Plane/Dynamo3DQG/Bx.hpp )
#include MAKE_STR( QUICC_MODEL_PATH/Boussinesq/Plane/Dynamo3DQG/By.hpp )
#include MAKE_STR( QUICC_MODEL_PATH/Boussinesq/Plane/Dynamo3DQG/fbx.hpp )
#include MAKE_STR( QUICC_MODEL_PATH/Boussinesq/Plane/Dynamo3DQG/fby.hpp )
#include MAKE_STR( QUICC_MODEL_PATH/Boussinesq/Plane/Dynamo3DQG/fbz.hpp )
#include "Enums/FieldIds.hpp"
#include "IoVariable/StateFileReader.hpp"
#include "IoVariable/StateFileWriter.hpp"
#include "IoVariable/VisualizationFileWriter.hpp"
#include "IoTools/IdToHuman.hpp"
#include "IoVariable/Cartesian1DNusseltZWriter.hpp"
#include "IoVariable/Cartesian1DScalarEnergyWriter.hpp"
#include "IoVariable/Cartesian1DStreamEnergyWriter.hpp"
#include "IoVariable/Cartesian1DMagneticEnergyWriter.hpp"
#include "IoStats/Cartesian1DScalarAvgWriter.hpp"
#include "Generator/States/RandomScalarState.hpp"
#include "Generator/States/CartesianExactScalarState.hpp"
#include "Generator/Visualizers/ScalarFieldVisualizer.hpp"
#include "Model/PhysicalModelBase.hpp"

namespace QuICC {

namespace Model {

namespace Boussinesq {

namespace Plane {

namespace Dynamo3DQG {

   const std::string PhysicalModel::PYMODULE = "boussinesq_dynamo3dqg";

   const std::string PhysicalModel::PYCLASS = "BoussinesqDynamo3DQG";

   void PhysicalModel::addEquations(SharedSimulation spSim)
   {
      // Add upright streamfunction equation
      spSim->addScalarEquation<Equations::Boussinesq::Plane::Dynamo3DQG::Streamfunction>();

      // Add upright vertical velocity equation
      spSim->addScalarEquation<Equations::Boussinesq::Plane::Dynamo3DQG::VelocityZ>();

      // Add upright transport equation
      spSim->addScalarEquation<Equations::Boussinesq::Plane::Dynamo3DQG::Transport>();


      // Add vertical vorticity equation
      spSim->addScalarEquation<Equations::Boussinesq::Plane::Dynamo3DQG::VorticityZ>();

      // Add mean heat computation
      spSim->addScalarEquation<Equations::Boussinesq::Plane::Dynamo3DQG::MeanHeat>();

      // Add Emfx computation
      spSim->addScalarEquation<Equations::Boussinesq::Plane::Dynamo3DQG::Emfx>();

      // Add Emfy computation
      spSim->addScalarEquation<Equations::Boussinesq::Plane::Dynamo3DQG::Emfy>();

      // Add Emfy computation
      spSim->addScalarEquation<Equations::Boussinesq::Plane::Dynamo3DQG::Bx>();

      // Add Emfy computation
      spSim->addScalarEquation<Equations::Boussinesq::Plane::Dynamo3DQG::By>();

      // Add fbx computation
      spSim->addScalarEquation<Equations::Boussinesq::Plane::Dynamo3DQG::fbx>();

      // Add fby computation
      spSim->addScalarEquation<Equations::Boussinesq::Plane::Dynamo3DQG::fby>();

      // Add fbz heat computation
      spSim->addScalarEquation<Equations::Boussinesq::Plane::Dynamo3DQG::fbz>();

      // Add fbz heat computation
      spSim->addScalarEquation<Equations::Boussinesq::Plane::Dynamo3DQG::Pressure>();
   }

   void PhysicalModel::addStates(SharedStateGenerator spGen)
   {
      // Generate "exact" solutions (trigonometric or monomial)
      if(false)
      {
         // Shared pointer to equation
         Equations::SharedCartesianExactScalarState spExact;

         // Add transport initial state generation equation
         spExact = spGen->addScalarEquation<Equations::CartesianExactScalarState>();
         spExact->setIdentity(PhysicalNames::TEMPERATURE);
         spExact->setStateType(Equations::CartesianExactStateIds::POLYCOSCOS);
         spExact->setModeOptions(1e0, 0.0, 1e0, 0.0, 1e0, 0.0);

         // Add streamfunction initial state generation equation
         spExact = spGen->addScalarEquation<Equations::CartesianExactScalarState>();
         spExact->setIdentity(PhysicalNames::STREAMFUNCTION);
         spExact->setStateType(Equations::CartesianExactStateIds::POLYCOSCOS);
         spExact->setModeOptions(2e0, 0.0, 1e0, 0.0, 1e0, 0.0);

         // Add vertical velocity initial state generation equation
         spExact = spGen->addScalarEquation<Equations::CartesianExactScalarState>();
         spExact->setIdentity(PhysicalNames::VELOCITYZ);
         spExact->setStateType(Equations::CartesianExactStateIds::POLYCOSCOS);
         spExact->setModeOptions(3e0, 0.0, 1e0, 0.0, 1e0, 0.0);

         // Add vertical velocity initial state generation equation
         spExact = spGen->addScalarEquation<Equations::CartesianExactScalarState>();
         spExact->setIdentity(PhysicalNames::DZ_MEANTEMPERATURE);
         spExact->setStateType(Equations::CartesianExactStateIds::POLYCOSCOS);
         spExact->setModeOptions(-1e0, 0.0, 1e0, 0.0, 1e0, 0.0);

         // Add BX initial state generation equation
         spExact = spGen->addScalarEquation<Equations::CartesianExactScalarState>();
         spExact->setIdentity(PhysicalNames::BX);
         spExact->setStateType(Equations::CartesianExactStateIds::POLYCOSCOS);
         spExact->setModeOptions(-1e0, 0.0, 1e0, 0.0, 1e0, 0.0);

         // Add BY initial state generation equation
         spExact = spGen->addScalarEquation<Equations::CartesianExactScalarState>();
         spExact->setIdentity(PhysicalNames::BY);
         spExact->setStateType(Equations::CartesianExactStateIds::POLYCOSCOS);
         spExact->setModeOptions(-1e0, 0.0, 1e0, 0.0, 1e0, 0.0);

         // Add BY initial state generation equation
         spExact = spGen->addScalarEquation<Equations::CartesianExactScalarState>();
         spExact->setIdentity(PhysicalNames::EMFY);
         spExact->setStateType(Equations::CartesianExactStateIds::POLYCOSCOS);
         spExact->setModeOptions(0, 0.0, 0, 0.0, 0, 0.0);
         
         // Add BY initial state generation equation
         spExact = spGen->addScalarEquation<Equations::CartesianExactScalarState>();
         spExact->setIdentity(PhysicalNames::EMFX);
         spExact->setStateType(Equations::CartesianExactStateIds::POLYCOSCOS);
         spExact->setModeOptions(0, 0.0, 0, 0.0, 0, 0.0);

         // Generate random spectrum
      } else
      {
         // Shared pointer to random initial state equation
         Equations::SharedRandomScalarState spRand;
         Equations::SharedCartesianExactScalarState spExact;

         // Add transport initial state generation equation
         spRand = spGen->addScalarEquation<Equations::RandomScalarState>();
         spRand->setIdentity(PhysicalNames::TEMPERATURE);
         spRand->setSpectrum(-1e-2, 1e-2, 1e4, 1e4, 1e4);

         // Add streamfunction initial state generation equation
         spRand = spGen->addScalarEquation<Equations::RandomScalarState>();
         spRand->setIdentity(PhysicalNames::STREAMFUNCTION);
         spRand->setSpectrum(-1e-2, 1e-2, 1e4, 1e4, 1e4);

         // Add vertical velocity initial state generation equation
         spRand = spGen->addScalarEquation<Equations::RandomScalarState>();
         spRand->setIdentity(PhysicalNames::VELOCITYZ);
         spRand->setSpectrum(-1e-2, 1e-2, 1e4, 1e4, 1e4);

         // Add vertical velocity initial state generation equation
         spExact = spGen->addScalarEquation<Equations::CartesianExactScalarState>();
         spExact->setIdentity(PhysicalNames::DZ_MEANTEMPERATURE);
         spExact->setStateType(Equations::CartesianExactStateIds::POLYCOSCOS);
         spExact->setModeOptions(-1e0, 0.0, 1e0, 0.0, 1e0, 0.0);

         // Add BX initial state generation equation
         spRand = spGen->addScalarEquation<Equations::RandomScalarState>();
         spRand->setIdentity(PhysicalNames::BX);
         spRand->setSpectrum(-1e-2, 1e-2, 1e4, 1e4, 1e4);

         // Add BY initial state generation equation
         spRand = spGen->addScalarEquation<Equations::RandomScalarState>();
         spRand->setIdentity(PhysicalNames::BY);
         spRand->setSpectrum(-1e-2, 1e-2, 1e4, 1e4, 1e4);

         // Add BY initial state generation equation
         spRand = spGen->addScalarEquation<Equations::RandomScalarState>();
         spRand->setIdentity(PhysicalNames::EMFY);
         spRand->setSpectrum(-1e-6, 1e-6, 1e4, 1e4, 1e4, Equations::RandomScalarState::ONLYMEAN);

         // Add BY initial state generation equation
         spRand = spGen->addScalarEquation<Equations::RandomScalarState>();
         spRand->setIdentity(PhysicalNames::EMFX);
         spRand->setSpectrum(-1e-6, 1e-6, 1e4, 1e4, 1e4,Equations::RandomScalarState::ONLYMEAN);

         // Add BY initial state generation equation
         spRand = spGen->addScalarEquation<Equations::RandomScalarState>();
         spRand->setIdentity(PhysicalNames::PRESSURE);
         spRand->setSpectrum(-1e-6, 1e-6, 1e4, 1e4, 1e4,Equations::RandomScalarState::ONLYMEAN);


      }

      // Add output file
      IoVariable::SharedStateFileWriter spOut(new IoVariable::StateFileWriter(SchemeType::type(), SchemeType::isRegular()));
      spOut->expect(PhysicalNames::TEMPERATURE);
      spOut->expect(PhysicalNames::STREAMFUNCTION);
      spOut->expect(PhysicalNames::VELOCITYZ);
      spOut->expect(PhysicalNames::DZ_MEANTEMPERATURE);
      spOut->expect(PhysicalNames::BX);
      spOut->expect(PhysicalNames::BY);
      spOut->expect(PhysicalNames::EMFX);
      spOut->expect(PhysicalNames::EMFY);
      spOut->expect(PhysicalNames::PRESSURE);
      spGen->addHdf5OutputFile(spOut);
   }

   void PhysicalModel::addVisualizers(SharedVisualizationGenerator spVis)
   {
      // Shared pointer to basic field visualizer
      Equations::SharedScalarFieldVisualizer spField;

      // Add transport field visualization
      spField = spVis->addScalarEquation<Equations::ScalarFieldVisualizer>();
      spField->setFields(true, false);
      spField->setIdentity(PhysicalNames::TEMPERATURE);

      // Add streamfunction field visualization
      spField = spVis->addScalarEquation<Equations::ScalarFieldVisualizer>();
      spField->setFields(true, false);
      spField->setIdentity(PhysicalNames::STREAMFUNCTION);

      // Add vertical velocity field visualization
      spField = spVis->addScalarEquation<Equations::ScalarFieldVisualizer>();
      spField->setFields(true, false);
      spField->setIdentity(PhysicalNames::VELOCITYZ);

      // Add background temperature profile visualization
      spField = spVis->addScalarEquation<Equations::ScalarFieldVisualizer>();
      spField->setFields(true, false);
      spField->setIdentity(PhysicalNames::DZ_MEANTEMPERATURE);

      // Add background temperature profile visualization
      spField = spVis->addScalarEquation<Equations::ScalarFieldVisualizer>();
      spField->setFields(true, false);
      spField->setIdentity(PhysicalNames::BX);

      // Add background temperature profile visualization
      spField = spVis->addScalarEquation<Equations::ScalarFieldVisualizer>();
      spField->setFields(true, false);
      spField->setIdentity(PhysicalNames::BY);

      // Add background temperature profile visualization
      spField = spVis->addScalarEquation<Equations::ScalarFieldVisualizer>();
      spField->setFields(true, false);
      spField->setIdentity(PhysicalNames::VORTICITYZ);

      // Add background temperature profile visualization
      spField = spVis->addScalarEquation<Equations::ScalarFieldVisualizer>();
      spField->setFields(true, false);
      spField->setIdentity(PhysicalNames::EMFX);

      // Add background temperature profile visualization
      spField = spVis->addScalarEquation<Equations::ScalarFieldVisualizer>();
      spField->setFields(true, false);
      spField->setIdentity(PhysicalNames::EMFY);

      // Add background temperature profile visualization
      spField = spVis->addScalarEquation<Equations::ScalarFieldVisualizer>();
      spField->setFields(true, false);
      spField->setIdentity(PhysicalNames::PRESSURE);


      // Add output file
      IoVariable::SharedVisualizationFileWriter spOut(new IoVariable::VisualizationFileWriter(SchemeType::type()));
      spOut->expect(PhysicalNames::TEMPERATURE);
      spOut->expect(PhysicalNames::STREAMFUNCTION);
      spOut->expect(PhysicalNames::VELOCITYZ);
      spOut->expect(PhysicalNames::DZ_MEANTEMPERATURE);
      spOut->expect(PhysicalNames::BX);
      spOut->expect(PhysicalNames::BY);
      spOut->expect(PhysicalNames::VORTICITYZ);
      spOut->expect(PhysicalNames::EMFY);
      spOut->expect(PhysicalNames::EMFX);
      spOut->expect(PhysicalNames::PRESSURE);
      spVis->addHdf5OutputFile(spOut);
   }

   void PhysicalModel::setVisualizationState(SharedVisualizationGenerator spVis)
   {
      // Create and add initial state file to IO
      IoVariable::SharedStateFileReader spIn(new IoVariable::StateFileReader("4Visu", SchemeType::type(), SchemeType::isRegular()));

      // Set expected fields
      spIn->expect(PhysicalNames::TEMPERATURE);
      spIn->expect(PhysicalNames::STREAMFUNCTION);
      spIn->expect(PhysicalNames::VELOCITYZ);
      spIn->expect(PhysicalNames::DZ_MEANTEMPERATURE);
      spIn->expect(PhysicalNames::BX);
      spIn->expect(PhysicalNames::BY);
      spIn->expect(PhysicalNames::VORTICITYZ);
      spIn->expect(PhysicalNames::EMFY);
      spIn->expect(PhysicalNames::EMFX);
      spIn->expect(PhysicalNames::PRESSURE);

      // Set simulation state
      spVis->setInitialState(spIn);
   }

   void PhysicalModel::addAsciiOutputFiles(SharedSimulation spSim)
   {
      // Create Nusselt number writer
      //IoVariable::SharedNusseltWriter spNusselt(new IoVariable::NusseltWriter(SchemeType::type()));
      //spNusselt->expect(PhysicalNames::DZ_MEANTEMPERATURE);
      //spSim->addAsciiOutputFile(spNusselt);

      // Create Nusselt number writer
      IoVariable::SharedCartesian1DNusseltZWriter spNusselt(new IoVariable::Cartesian1DNusseltZWriter(SchemeType::type()));
      spNusselt->expect(PhysicalNames::DZ_MEANTEMPERATURE);
      spSim->addAsciiOutputFile(spNusselt);

      // Create temperature energy writer
      IoVariable::SharedCartesian1DScalarEnergyWriter spTemp(new IoVariable::Cartesian1DScalarEnergyWriter("temperature", SchemeType::type()));
      spTemp->expect(PhysicalNames::TEMPERATURE);
      spSim->addAsciiOutputFile(spTemp);

      // Create kinetic energy writer
      IoVariable::SharedCartesian1DStreamEnergyWriter spStream(new IoVariable::Cartesian1DStreamEnergyWriter("kinetic", SchemeType::type()));
      spStream->expect(PhysicalNames::STREAMFUNCTION);
      spStream->expect(PhysicalNames::VELOCITYZ);
      spSim->addAsciiOutputFile(spStream);

      // Create magnetic energy writer
      IoVariable::SharedCartesian1DMagneticEnergyWriter spMag(new IoVariable::Cartesian1DMagneticEnergyWriter("magnetic", SchemeType::type()));
      spMag->expect(PhysicalNames::BX);
      spMag->expect(PhysicalNames::BY);
      spSim->addAsciiOutputFile(spMag);

      // Create temperature energy writer
      //   IoVariable::SharedCartesian1DScalarEnergyWriter spMag(new IoVariable::Cartesian1DScalarEnergyWriter("magX", SchemeType::type()));
      //   spMag->expect(PhysicalNames::BX);
      //   spSim->addAsciiOutputFile(spMag);

   }

   void PhysicalModel::addStatsOutputFiles(SharedSimulation spSim)
   {
      // Create Avg temperature writer
      IoStats::SharedCartesian1DScalarAvgWriter spAvg(new IoStats::Cartesian1DScalarAvgWriter("temperature",SchemeType::type()));
      spAvg->expect(PhysicalNames::TEMPERATURE);
      spSim->addStatsOutputFile(spAvg);

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

      // Add mean temperature to ouput file
      spState->expect(PhysicalNames::DZ_MEANTEMPERATURE);
      spState->expect(PhysicalNames::VORTICITYZ);
      spState->expect(PhysicalNames::EMFY);
      spState->expect(PhysicalNames::EMFX);
      spState->expect(PhysicalNames::PRESSURE);

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

      spInit->expect(PhysicalNames::EMFY);
      spInit->expect(PhysicalNames::EMFX);
      spInit->expect(PhysicalNames::PRESSURE);

      // Set simulation state
      spSim->setInitialState(spInit);
   }

}
}
}
}
}
