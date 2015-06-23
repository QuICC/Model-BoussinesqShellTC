/** 
 * @file BoussinesqNoTiltedFPlane3DQGModel.cpp
 * @brief Source of the Boussinesq tilted F-plane 3DQG physical model in nonorthogonal coordinates
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
#include "PhysicalModels/BoussinesqNoTiltedFPlane3DQGModel.hpp"

// Project includes
//
#include "Enums/FieldIds.hpp"
#include "IoVariable/StateFileReader.hpp"
#include "IoVariable/StateFileWriter.hpp"
#include "IoVariable/VisualizationFileWriter.hpp"
#include "IoTools/IdToHuman.hpp"
#include "Equations/Asymptotics/FPlane3DQG/Boussinesq/BoussinesqNoTiltedFPlane3DQGTransport.hpp"
#include "Equations/Asymptotics/FPlane3DQG/Boussinesq/BoussinesqNoTiltedFPlane3DQGNoStreamfunction.hpp"
#include "Equations/Asymptotics/FPlane3DQG/Boussinesq/BoussinesqNoTiltedFPlane3DQGNoVelocityZ.hpp"
#include "Equations/Asymptotics/FPlane3DQG/Boussinesq/BoussinesqNoTiltedFPlane3DQGNoVorticityZ.hpp"
#include "Equations/Asymptotics/FPlane3DQG/Boussinesq/BoussinesqNoTiltedFPlane3DQGMeanHeat.hpp"
#include "IoVariable/NusseltWriter.hpp"
#include "IoVariable/Cartesian1DScalarEnergyWriter.hpp"
#include "IoVariable/Cartesian1DStreamEnergyWriter.hpp"
#include "Generator/States/RandomScalarState.hpp"
#include "Generator/States/CartesianExactScalarState.hpp"
#include "Generator/Visualizers/ScalarFieldVisualizer.hpp"
#include "Generator/Visualizers/TiltedScalarFieldVisualizer.hpp"
#include "PhysicalModels/PhysicalModelBase.hpp"

namespace GeoMHDiSCC {

   const std::string BoussinesqNoTiltedFPlane3DQGModel::PYMODULE = "boussinesq_notilted_fplane3dqg";

   const std::string BoussinesqNoTiltedFPlane3DQGModel::PYCLASS = "BoussinesqNoTiltedFPlane3DQG";

   void BoussinesqNoTiltedFPlane3DQGModel::addEquations(SharedSimulation spSim)
   {
      // Add non orthogonal streamfunction equation
      spSim->addScalarEquation<Equations::BoussinesqNoTiltedFPlane3DQGNoStreamfunction>();
      
      // Add non orthogonal vertical velocity equation
      spSim->addScalarEquation<Equations::BoussinesqNoTiltedFPlane3DQGNoVelocityZ>();
      
      // Add upright transport equation
      spSim->addScalarEquation<Equations::BoussinesqNoTiltedFPlane3DQGTransport>();

      // Add non orthogonal vertical vorticity equation
      spSim->addScalarEquation<Equations::BoussinesqNoTiltedFPlane3DQGNoVorticityZ>();

      
      // Add mean heat computation
      spSim->addScalarEquation<Equations::BoussinesqNoTiltedFPlane3DQGMeanHeat>();
   }

   void BoussinesqNoTiltedFPlane3DQGModel::addStates(SharedStateGenerator spGen)
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
         spExact->setIdentity(PhysicalNames::NO_STREAMFUNCTION);
         spExact->setStateType(Equations::CartesianExactStateIds::POLYCOSCOS);
         spExact->setModeOptions(2e0, 3.0, 1e0, 3.0, 1e0, 3.0);

         // Add vertical velocity initial state generation equation
         spExact = spGen->addScalarEquation<Equations::CartesianExactScalarState>();
         spExact->setIdentity(PhysicalNames::NO_VELOCITYZ);
         spExact->setStateType(Equations::CartesianExactStateIds::POLYCOSCOS);
         spExact->setModeOptions(1e0, -2.0, 1e0, 5.0, 7e0, 4.0);

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
         spRand->setIdentity(PhysicalNames::NO_STREAMFUNCTION);
         spRand->setSpectrum(-1.0e-4, 1.0e-4, 1e4, 1e4, 1e4);

         // Add vertical velocity initial state generation equation
         spRand = spGen->addScalarEquation<Equations::RandomScalarState>();
         spRand->setIdentity(PhysicalNames::NO_VELOCITYZ);
         spRand->setSpectrum(-1.0e-4, 1.0e-4, 1e4, 1e4, 1e4);
      }

      // Add output file
      IoVariable::SharedStateFileWriter spOut(new IoVariable::StateFileWriter(SchemeType::type(), SchemeType::isRegular()));
      spOut->expect(PhysicalNames::TEMPERATURE);
      spOut->expect(PhysicalNames::NO_STREAMFUNCTION);
      spOut->expect(PhysicalNames::NO_VELOCITYZ);
      spGen->addHdf5OutputFile(spOut);
   }

   void BoussinesqNoTiltedFPlane3DQGModel::addVisualizers(SharedVisualizationGenerator spVis)
   {
      // Shared pointer to basic field visualizer
      Equations::SharedScalarFieldVisualizer spField;

      // Add transport field visualization
      spField = spVis->addScalarEquation<Equations::ScalarFieldVisualizer>();
      spField->setFields(true, false);
      spField->setIdentity(PhysicalNames::TEMPERATURE);
      
      // Add non orthogonal streamfunction field visualization
      spField = spVis->addScalarEquation<Equations::ScalarFieldVisualizer>();
      spField->setFields(true, false);
      spField->setIdentity(PhysicalNames::NO_STREAMFUNCTION);
      
      // Add non orthogonal vertical velocity field visualization
      spField = spVis->addScalarEquation<Equations::ScalarFieldVisualizer>();
      spField->setFields(true, false);
      spField->setIdentity(PhysicalNames::NO_VELOCITYZ);
      
      // Add mean temperature z gradient field visualization
      spField = spVis->addScalarEquation<Equations::ScalarFieldVisualizer>();
      spField->setFields(true, false);
      spField->setIdentity(PhysicalNames::DZ_MEANTEMPERATURE);
      
      // Add vertical vorticity field visualization
      spField = spVis->addScalarEquation<Equations::ScalarFieldVisualizer>();
      spField->setFields(true, false);
      spField->setIdentity(PhysicalNames::NO_VORTICITYZ);

      // Add output file
      IoVariable::SharedVisualizationFileWriter spOut(new IoVariable::VisualizationFileWriter(SchemeType::type()));
      spOut->expect(PhysicalNames::TEMPERATURE);
      spOut->expect(PhysicalNames::NO_STREAMFUNCTION);
      spOut->expect(PhysicalNames::NO_VELOCITYZ);
      spOut->expect(PhysicalNames::DZ_MEANTEMPERATURE);
      spOut->expect(PhysicalNames::NO_VORTICITYZ);
      spVis->addHdf5OutputFile(spOut);
   }

   void BoussinesqNoTiltedFPlane3DQGModel::setVisualizationState(SharedVisualizationGenerator spVis)
   {
      // Create and add initial state file to IO
      IoVariable::SharedStateFileReader spIn(new IoVariable::StateFileReader("4Visu", SchemeType::type(), SchemeType::isRegular()));

      // Set expected fields
      spIn->expect(PhysicalNames::TEMPERATURE);
      spIn->expect(PhysicalNames::NO_STREAMFUNCTION);
      spIn->expect(PhysicalNames::NO_VELOCITYZ);
      spIn->expect(PhysicalNames::DZ_MEANTEMPERATURE);

      // Set simulation state
      spVis->setInitialState(spIn);
   }

   void BoussinesqNoTiltedFPlane3DQGModel::addAsciiOutputFiles(SharedSimulation spSim)
   {
      // Create Nusselt number writer
      IoVariable::SharedNusseltWriter spState(new IoVariable::NusseltWriter(SchemeType::type()));
      spState->expect(PhysicalNames::DZ_MEANTEMPERATURE);
      spSim->addAsciiOutputFile(spState);

      // Create temperature energy writer
      IoVariable::SharedCartesian1DScalarEnergyWriter spTemp(new IoVariable::Cartesian1DScalarEnergyWriter("temperature", SchemeType::type()));
      spTemp->expect(PhysicalNames::TEMPERATURE);
      spSim->addAsciiOutputFile(spTemp);

      // Create kinetic energy writer
      IoVariable::SharedCartesian1DStreamEnergyWriter spStream(new IoVariable::Cartesian1DStreamEnergyWriter("kinetic", SchemeType::type()));
      spStream->expect(PhysicalNames::NO_STREAMFUNCTION);
      spStream->expect(PhysicalNames::NO_VELOCITYZ);
      spSim->addAsciiOutputFile(spStream);
   }

   void BoussinesqNoTiltedFPlane3DQGModel::addHdf5OutputFiles(SharedSimulation spSim)
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

   void BoussinesqNoTiltedFPlane3DQGModel::setInitialState(SharedSimulation spSim)
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
