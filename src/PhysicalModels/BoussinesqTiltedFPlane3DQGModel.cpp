/** 
 * @file BoussinesqTiltedFPlane3DQGModel.cpp
 * @brief Source of the Boussinesq tilted F-plane 3DQG physical model
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
#include "PhysicalModels/BoussinesqTiltedFPlane3DQGModel.hpp"

// Project includes
//
#include "Enums/FieldIds.hpp"
#include "IoVariable/StateFileReader.hpp"
#include "IoVariable/StateFileWriter.hpp"
#include "IoVariable/VisualizationFileWriter.hpp"
#include "IoTools/IdToHuman.hpp"
#include "Equations/Asymptotics/FPlane3DQG/Boussinesq/BoussinesqTiltedFPlane3DQGStreamfunction.hpp"
#include "Equations/Asymptotics/FPlane3DQG/Boussinesq/BoussinesqTiltedFPlane3DQGVelocityZ.hpp"
#include "Equations/Asymptotics/FPlane3DQG/Boussinesq/BoussinesqTiltedFPlane3DQGTransport.hpp"
#include "Equations/Asymptotics/FPlane3DQG/Boussinesq/BoussinesqTiltedFPlane3DQGNoStreamfunction.hpp"
#include "Equations/Asymptotics/FPlane3DQG/Boussinesq/BoussinesqTiltedFPlane3DQGNoVelocityZ.hpp"
#include "Equations/Asymptotics/FPlane3DQG/Boussinesq/BoussinesqTiltedFPlane3DQGNoVorticityZ.hpp"
#include "Equations/Asymptotics/FPlane3DQG/Boussinesq/BoussinesqTiltedFPlane3DQGMeanHeat.hpp"
#include "IoVariable/Cartesian1DNusseltZWriter.hpp"
#include "IoVariable/Cartesian1DScalarEnergyWriter.hpp"
#include "IoVariable/Cartesian1DStreamEnergyWriter.hpp"
#include "Generator/States/RandomScalarState.hpp"
#include "Generator/States/CartesianExactScalarState.hpp"
#include "Generator/Visualizers/ScalarFieldVisualizer.hpp"
#include "Generator/Visualizers/TiltedScalarFieldVisualizer.hpp"
#include "PhysicalModels/PhysicalModelBase.hpp"

namespace GeoMHDiSCC {

   const std::string BoussinesqTiltedFPlane3DQGModel::PYMODULE = "boussinesq_tilted_fplane3dqg_r";

   const std::string BoussinesqTiltedFPlane3DQGModel::PYCLASS = "BoussinesqTiltedFPlane3DQG";

   void BoussinesqTiltedFPlane3DQGModel::addEquations(SharedSimulation spSim)
   {
      // Add upright streamfunction equation
      spSim->addScalarEquation<Equations::BoussinesqTiltedFPlane3DQGStreamfunction>();
      
      // Add upright vertical velocity equation
      spSim->addScalarEquation<Equations::BoussinesqTiltedFPlane3DQGVelocityZ>();
      
      // Add upright transport equation
      spSim->addScalarEquation<Equations::BoussinesqTiltedFPlane3DQGTransport>();


      // Add non orthogonal streamfunction equation
      spSim->addScalarEquation<Equations::BoussinesqTiltedFPlane3DQGNoStreamfunction>(SolveTiming::BEFORE);

      // Add non orthogonal vertical velocity equation
      spSim->addScalarEquation<Equations::BoussinesqTiltedFPlane3DQGNoVelocityZ>(SolveTiming::BEFORE);

      // Add non orthogonal streamfunction equation
      spSim->addScalarEquation<Equations::BoussinesqTiltedFPlane3DQGNoStreamfunction>(SolveTiming::AFTER);

      // Add non orthogonal vertical velocity equation
      spSim->addScalarEquation<Equations::BoussinesqTiltedFPlane3DQGNoVelocityZ>(SolveTiming::AFTER);

      // Add non orthogonal vertical vorticity equation
      spSim->addScalarEquation<Equations::BoussinesqTiltedFPlane3DQGNoVorticityZ>();

      
      // Add mean heat computation
      spSim->addScalarEquation<Equations::BoussinesqTiltedFPlane3DQGMeanHeat>();
   }

   void BoussinesqTiltedFPlane3DQGModel::addStates(SharedStateGenerator spGen)
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
         spExact->setModeOptions(1e0, 4.0, -1e0, 4.0, 1e0, 4.0);

         // Add streamfunction initial state generation equation
         spExact = spGen->addScalarEquation<Equations::CartesianExactScalarState>();
         spExact->setIdentity(PhysicalNames::STREAMFUNCTION);
         spExact->setStateType(Equations::CartesianExactStateIds::POLYCOSCOS);
         spExact->setModeOptions(2e0, 3.0, 1e0, 3.0, 1e0, 3.0);

         // Add vertical velocity initial state generation equation
         spExact = spGen->addScalarEquation<Equations::CartesianExactScalarState>();
         spExact->setIdentity(PhysicalNames::VELOCITYZ);
         spExact->setStateType(Equations::CartesianExactStateIds::POLYCOSCOS);
         spExact->setModeOptions(1e0, -2.0, 1e0, 5.0, 7e0, 4.0);

         // Add initial mean temperature gradient state generation equation
         spExact = spGen->addScalarEquation<Equations::CartesianExactScalarState>();
         spExact->setIdentity(PhysicalNames::DZ_MEANTEMPERATURE);
         spExact->setStateType(Equations::CartesianExactStateIds::POLYCOSCOS);
         spExact->setModeOptions(-1e0, 0.0, 1e0, 0.0, 1e0, 0.0);

      // Generate random spectrum
      } else
      {
         // Shared pointer to random initial state equation
         Equations::SharedRandomScalarState spRand;

         // Add transport initial state generation equation
         spRand = spGen->addScalarEquation<Equations::RandomScalarState>();
         spRand->setIdentity(PhysicalNames::TEMPERATURE);
         spRand->setSpectrum(-1.0e-4, 1.0e-4, 1e4, 1e4, 1e4);

         // Add streamfunction initial state generation equation
         spRand = spGen->addScalarEquation<Equations::RandomScalarState>();
         spRand->setIdentity(PhysicalNames::STREAMFUNCTION);
         spRand->setSpectrum(-1.0e-4, 1.0e-4, 1e4, 1e4, 1e4);

         // Add vertical velocity initial state generation equation
         spRand = spGen->addScalarEquation<Equations::RandomScalarState>();
         spRand->setIdentity(PhysicalNames::VELOCITYZ);
         spRand->setSpectrum(-1.0e-4, 1.0e-4, 1e4, 1e4, 1e4);

         // Add initial mean temperature gradient state generation equation
         Equations::SharedCartesianExactScalarState spExact;
         spExact = spGen->addScalarEquation<Equations::CartesianExactScalarState>();
         spExact->setIdentity(PhysicalNames::DZ_MEANTEMPERATURE);
         spExact->setStateType(Equations::CartesianExactStateIds::POLYCOSCOS);
         spExact->setModeOptions(-1e0, 0.0, 1e0, 0.0, 1e0, 0.0);
      }

      // Add output file
      IoVariable::SharedStateFileWriter spOut(new IoVariable::StateFileWriter(SchemeType::type(), SchemeType::isRegular()));
      spOut->expect(PhysicalNames::TEMPERATURE);
      spOut->expect(PhysicalNames::STREAMFUNCTION);
      spOut->expect(PhysicalNames::VELOCITYZ);
      spOut->expect(PhysicalNames::DZ_MEANTEMPERATURE);
      spGen->addHdf5OutputFile(spOut);
   }

   void BoussinesqTiltedFPlane3DQGModel::addVisualizers(SharedVisualizationGenerator spVis)
   {
      // Shared pointer to basic field visualizer
      Equations::SharedScalarFieldVisualizer spField;
      Equations::SharedTiltedScalarFieldVisualizer spTilted;

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
      
      // Add vertical velocity field visualization
      spField = spVis->addScalarEquation<Equations::ScalarFieldVisualizer>();
      spField->setFields(true, false);
      spField->setIdentity(PhysicalNames::DZ_MEANTEMPERATURE);
//
//      // Add transport field visualization
//      spTilted = spVis->addScalarEquation<Equations::TiltedScalarFieldVisualizer>();
//      spTilted->setFields(true, false);
//      spTilted->setDataField(PhysicalNames::TEMPERATURE);
//      spTilted->setIdentity(PhysicalNames::TILTED_TEMPERATURE);
//      
//      // Add nonorthogonal streamfunction field visualization
//      spTilted = spVis->addScalarEquation<Equations::TiltedScalarFieldVisualizer>();
//      spTilted->setFields(true, false);
//      spTilted->setDataField(PhysicalNames::STREAMFUNCTION);
//      spTilted->setIdentity(PhysicalNames::TILTED_STREAMFUNCTION);
//      
//      // Add nonorthogonal vertical velocity field visualization
//      spTilted = spVis->addScalarEquation<Equations::TiltedScalarFieldVisualizer>();
//      spTilted->setFields(true, false);
//      spTilted->setDataField(PhysicalNames::VELOCITYZ);
//      spTilted->setIdentity(PhysicalNames::TILTED_VELOCITYZ);
//
//      // Add non orthogonal streamfunction equation
//      spVis->addScalarEquation<Equations::BoussinesqTiltedFPlane3DQGNoStreamfunction>(SolveTiming::AFTER);
//      // Add non orthogonal vertical velocity equation
//      spVis->addScalarEquation<Equations::BoussinesqTiltedFPlane3DQGNoVelocityZ>(SolveTiming::AFTER);
//      // Add non orthogonal vertical vorticity equation
//      spVis->addScalarEquation<Equations::BoussinesqTiltedFPlane3DQGNoVorticityZ>();
//      
//      // Add streamfunction field visualization
//      spField = spVis->addScalarEquation<Equations::ScalarFieldVisualizer>();
//      spField->setFields(true, false);
//      spField->setIdentity(PhysicalNames::NO_STREAMFUNCTION);
//      
//      // Add vertical velocity field visualization
//      spField = spVis->addScalarEquation<Equations::ScalarFieldVisualizer>();
//      spField->setFields(true, false);
//      spField->setIdentity(PhysicalNames::NO_VELOCITYZ);
//      
//      // Add nonorthogonal streamfunction field visualization
//      spTilted = spVis->addScalarEquation<Equations::TiltedScalarFieldVisualizer>();
//      spTilted->setFields(true, false);
//      spTilted->setDataField(PhysicalNames::NO_STREAMFUNCTION);
//      spTilted->setIdentity(PhysicalNames::TILTED_NO_STREAMFUNCTION);
//      
//      // Add nonorthogonal vertical velocity field visualization
//      spTilted = spVis->addScalarEquation<Equations::TiltedScalarFieldVisualizer>();
//      spTilted->setFields(true, false);
//      spTilted->setDataField(PhysicalNames::NO_VELOCITYZ);
//      spTilted->setIdentity(PhysicalNames::TILTED_NO_VELOCITYZ);
//      
//      // Add nonorthogonal vertical vorticity field visualization
//      spTilted = spVis->addScalarEquation<Equations::TiltedScalarFieldVisualizer>();
//      spTilted->setFields(true, false);
//      spTilted->setDataField(PhysicalNames::NO_VORTICITYZ);
//      spTilted->setIdentity(PhysicalNames::TILTED_NO_VORTICITYZ);

      // Add output file
      IoVariable::SharedVisualizationFileWriter spOut(new IoVariable::VisualizationFileWriter(SchemeType::type()));
      spOut->expect(PhysicalNames::TEMPERATURE);
      spOut->expect(PhysicalNames::STREAMFUNCTION);
      spOut->expect(PhysicalNames::VELOCITYZ);
      spOut->expect(PhysicalNames::DZ_MEANTEMPERATURE);
//      spOut->expect(PhysicalNames::TILTED_STREAMFUNCTION);
//      spOut->expect(PhysicalNames::TILTED_VELOCITYZ);
//      spOut->expect(PhysicalNames::NO_STREAMFUNCTION);
//      spOut->expect(PhysicalNames::NO_VELOCITYZ);
//      spOut->expect(PhysicalNames::NO_VORTICITYZ);
//      spOut->expect(PhysicalNames::TILTED_NO_STREAMFUNCTION);
//      spOut->expect(PhysicalNames::TILTED_NO_VELOCITYZ);
//      spOut->expect(PhysicalNames::TILTED_NO_VORTICITYZ);
      spVis->addHdf5OutputFile(spOut);
   }

   void BoussinesqTiltedFPlane3DQGModel::setVisualizationState(SharedVisualizationGenerator spVis)
   {
      // Create and add initial state file to IO
      IoVariable::SharedStateFileReader spIn(new IoVariable::StateFileReader("4Visu", SchemeType::type(), SchemeType::isRegular()));

      // Set expected fields
      spIn->expect(PhysicalNames::TEMPERATURE);
      spIn->expect(PhysicalNames::STREAMFUNCTION);
      spIn->expect(PhysicalNames::VELOCITYZ);
      spIn->expect(PhysicalNames::DZ_MEANTEMPERATURE);

      // Set simulation state
      spVis->setInitialState(spIn);
   }

   void BoussinesqTiltedFPlane3DQGModel::addAsciiOutputFiles(SharedSimulation spSim)
   {
      // Create Nusselt number writer
      IoVariable::SharedCartesian1DNusseltZWriter spState(new IoVariable::Cartesian1DNusseltZWriter(SchemeType::type()));
      spState->expect(PhysicalNames::DZ_MEANTEMPERATURE);
      spSim->addAsciiOutputFile(spState);

      // Create temperature energy writer
      IoVariable::SharedCartesian1DScalarEnergyWriter spTemp(new IoVariable::Cartesian1DScalarEnergyWriter("temperature", SchemeType::type()));
      spTemp->expect(PhysicalNames::TEMPERATURE);
      spSim->addAsciiOutputFile(spTemp);

      // Create kinetic energy writer
      IoVariable::SharedCartesian1DStreamEnergyWriter spStream(new IoVariable::Cartesian1DStreamEnergyWriter("kinetic", SchemeType::type()));
      spStream->expect(PhysicalNames::STREAMFUNCTION);
      spStream->expect(PhysicalNames::VELOCITYZ);
      spSim->addAsciiOutputFile(spStream);
   }

   void BoussinesqTiltedFPlane3DQGModel::addHdf5OutputFiles(SharedSimulation spSim)
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

      spSim->addHdf5OutputFile(spState);
   }

   void BoussinesqTiltedFPlane3DQGModel::addStatsOutputFiles(SharedSimulation spSim)
   {
   }

   void BoussinesqTiltedFPlane3DQGModel::setInitialState(SharedSimulation spSim)
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
