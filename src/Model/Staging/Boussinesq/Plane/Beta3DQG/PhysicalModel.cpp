/** 
 * @file BoussinesqBeta3DQGPerModel.cpp
 * @brief Source of the Boussinesq F-plane 3DQG physical model
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
#include "PhysicalModels/BoussinesqBeta3DQGPerModel.hpp"

// Project includes
//
#include "Enums/FieldIds.hpp"
#include "IoVariable/StateFileReader.hpp"
#include "IoVariable/StateFileWriter.hpp"
#include "IoVariable/VisualizationFileWriter.hpp"
#include "IoTools/IdToHuman.hpp"
#include "Equations/Asymptotics/Beta3DQG/Boussinesq/BoussinesqBeta3DQGPerStreamfunction.hpp"
#include "Equations/Asymptotics/Beta3DQG/Boussinesq/BoussinesqBeta3DQGPerVelocityZ.hpp"
#include "Equations/Asymptotics/Beta3DQG/Boussinesq/BoussinesqBeta3DQGPerTransport.hpp"
#include "Equations/Asymptotics/Beta3DQG/Boussinesq/BoussinesqBeta3DQGPerVorticityZ.hpp"
#include "Equations/Asymptotics/Beta3DQG/Boussinesq/BoussinesqBeta3DQGPerMeanHeat.hpp"
#include "IoVariable/Cartesian1DNusseltXWriter.hpp"
#include "IoVariable/Cartesian1DScalarEnergyWriter.hpp"
#include "IoVariable/Cartesian1DStreamEnergyWriter.hpp"
#include "Generator/States/RandomScalarState.hpp"
#include "Generator/States/CartesianExactScalarState.hpp"
#include "Generator/Visualizers/ScalarFieldVisualizer.hpp"
#include "PhysicalModels/PhysicalModelBase.hpp"

namespace QuICC {

   const std::string BoussinesqBeta3DQGPerModel::PYMODULE = "boussinesq_beta3dqg_per";

   const std::string BoussinesqBeta3DQGPerModel::PYCLASS = "BoussinesqBeta3DQGPer";

   void BoussinesqBeta3DQGPerModel::addEquations(SharedSimulation spSim)
   {
      // Add streamfunction equation
      spSim->addScalarEquation<Equations::BoussinesqBeta3DQGPerStreamfunction>();
      
      // Add vertical velocity computation
      spSim->addScalarEquation<Equations::BoussinesqBeta3DQGPerVelocityZ>();
      
      // Add transport equation
      spSim->addScalarEquation<Equations::BoussinesqBeta3DQGPerTransport>();
      
      // Add vertical vorticity computation
      spSim->addScalarEquation<Equations::BoussinesqBeta3DQGPerVorticityZ>();

      
      // Add mean heat computation
      spSim->addScalarEquation<Equations::BoussinesqBeta3DQGPerMeanHeat>();
   }

   void BoussinesqBeta3DQGPerModel::addStates(SharedStateGenerator spGen)
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
         spExact->setModeOptions(1e0, 3.0, 1e0, 3.0, 1e0, 3.0);

         // Add streamfunction initial state generation equation
         spExact = spGen->addScalarEquation<Equations::CartesianExactScalarState>();
         spExact->setIdentity(PhysicalNames::STREAMFUNCTION);
         spExact->setStateType(Equations::CartesianExactStateIds::POLYCOSCOS);
         spExact->setModeOptions(1e0, 7.0, 1e0, 5.0, 1e0, 4.0);

         // Add vertical velocity initial state generation equation
         spExact = spGen->addScalarEquation<Equations::CartesianExactScalarState>();
         spExact->setIdentity(PhysicalNames::VELOCITYZ);
         spExact->setStateType(Equations::CartesianExactStateIds::POLYCOSCOS);
         spExact->setModeOptions(1e0, 5.0, 1e0, 4.0, 1e0, 7.0);

         // Add vertical vorticity initial state generation equation
         spExact = spGen->addScalarEquation<Equations::CartesianExactScalarState>();
         spExact->setIdentity(PhysicalNames::VORTICITYZ);
         spExact->setStateType(Equations::CartesianExactStateIds::POLYCOSCOS);
         spExact->setModeOptions(1e0, 1.0, 1e0, 1.0, 1e0, 1.0);

      // Generate random spectrum
      } else
      {
         // Shared pointer to random initial state equation
         Equations::SharedRandomScalarState spRand;

         // Add transport initial state generation equation
         spRand = spGen->addScalarEquation<Equations::RandomScalarState>();
         spRand->setIdentity(PhysicalNames::TEMPERATURE);
         spRand->setSpectrum(-1e-6, 1e-6, 1e4, 1e4, 1e4);

         // Add streamfunction initial state generation equation
         spRand = spGen->addScalarEquation<Equations::RandomScalarState>();
         spRand->setIdentity(PhysicalNames::STREAMFUNCTION);
         spRand->setSpectrum(-1e-6, 1e-6, 1e4, 1e4, 1e4);

         // Add vertical velocity initial state generation equation
         spRand = spGen->addScalarEquation<Equations::RandomScalarState>();
         spRand->setIdentity(PhysicalNames::VELOCITYZ);
         spRand->setSpectrum(-1e-6, 1e-6, 1e4, 1e4, 1e4);

         // Add vertical vorticity initial state generation equation
         spRand = spGen->addScalarEquation<Equations::RandomScalarState>();
         spRand->setIdentity(PhysicalNames::VORTICITYZ);
         spRand->setSpectrum(-1e-6, 1e-6, 1e4, 1e4, 1e4);
      }

      // Add output file
      IoVariable::SharedStateFileWriter spOut(new IoVariable::StateFileWriter(SchemeType::type(), SchemeType::isRegular()));
      spOut->expect(PhysicalNames::TEMPERATURE);
      spOut->expect(PhysicalNames::STREAMFUNCTION);
      spOut->expect(PhysicalNames::VELOCITYZ);
      spOut->expect(PhysicalNames::VORTICITYZ);
      spGen->addHdf5OutputFile(spOut);
   }

   void BoussinesqBeta3DQGPerModel::addVisualizers(SharedVisualizationGenerator spVis)
   {
      // Shared pointer to basic field visualizer
      Equations::SharedScalarFieldVisualizer spField;

      // Add transport field visualization
      spField = spVis->addScalarEquation<Equations::ScalarFieldVisualizer>();
      spField->setFields(true, false);
      spField->setIdentity(PhysicalNames::TEMPERATURE);
      
      // Add streamfunction field visualization
      spField = spVis->addScalarEquation<Equations::ScalarFieldVisualizer>();
      spField->setFields(true, true);
      spField->setIdentity(PhysicalNames::STREAMFUNCTION);
      
      // Add vertical velocity field visualization
      spField = spVis->addScalarEquation<Equations::ScalarFieldVisualizer>();
      spField->setFields(true, false);
      spField->setIdentity(PhysicalNames::VELOCITYZ);
      
//      // Add vertical velocity field visualization
//      spField = spVis->addScalarEquation<Equations::ScalarFieldVisualizer>();
//      spField->setFields(true, false);
//      spField->setIdentity(PhysicalNames::VORTICITYZ);
      
      // Add vertical velocity field visualization
      spField = spVis->addScalarEquation<Equations::ScalarFieldVisualizer>();
      spField->setFields(true, false);
      spField->setIdentity(PhysicalNames::DX_MEANTEMPERATURE);

      // Add output file
      IoVariable::SharedVisualizationFileWriter spOut(new IoVariable::VisualizationFileWriter(SchemeType::type()));
      spOut->expect(PhysicalNames::TEMPERATURE);
      spOut->expect(PhysicalNames::STREAMFUNCTION);
      spOut->expect(PhysicalNames::VELOCITYZ);
//      spOut->expect(PhysicalNames::VORTICITYZ);
      spOut->expect(PhysicalNames::DX_MEANTEMPERATURE);
      spVis->addHdf5OutputFile(spOut);
   }

   void BoussinesqBeta3DQGPerModel::setVisualizationState(SharedVisualizationGenerator spVis)
   {
      // Create and add initial state file to IO
      IoVariable::SharedStateFileReader spIn(new IoVariable::StateFileReader("4Visu", SchemeType::type(), SchemeType::isRegular()));

      // Set expected fields
      spIn->expect(PhysicalNames::TEMPERATURE);
      spIn->expect(PhysicalNames::STREAMFUNCTION);
      spIn->expect(PhysicalNames::VELOCITYZ);
      spIn->expect(PhysicalNames::VORTICITYZ);
      spIn->expect(PhysicalNames::DX_MEANTEMPERATURE);

      // Set simulation state
      spVis->setInitialState(spIn);
   }

   void BoussinesqBeta3DQGPerModel::addAsciiOutputFiles(SharedSimulation spSim)
   {
      // Create Nusselt number writer
      IoVariable::SharedCartesian1DNusseltXWriter spNusselt(new IoVariable::Cartesian1DNusseltXWriter(SchemeType::type()));
      spNusselt->expect(PhysicalNames::TEMPERATURE);
      spNusselt->expect(PhysicalNames::DX_MEANTEMPERATURE);
      spSim->addAsciiOutputFile(spNusselt);

      // Create temperature energy writer
      IoVariable::SharedCartesian1DScalarEnergyWriter spTemp(new IoVariable::Cartesian1DScalarEnergyWriter("temperature", SchemeType::type()));
      spTemp->expect(PhysicalNames::TEMPERATURE);
      spSim->addAsciiOutputFile(spTemp);

      // Create kinetic energy writer
      IoVariable::SharedCartesian1DStreamEnergyWriter spStream(new IoVariable::Cartesian1DStreamEnergyWriter("kinetic", SchemeType::type(), false, true));
      spStream->expect(PhysicalNames::STREAMFUNCTION);
      spStream->expect(PhysicalNames::VELOCITYZ);
      spSim->addAsciiOutputFile(spStream);
   }

   void BoussinesqBeta3DQGPerModel::addHdf5OutputFiles(SharedSimulation spSim)
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
      spState->expect(PhysicalNames::DX_MEANTEMPERATURE);

      spSim->addHdf5OutputFile(spState);
   }

   void BoussinesqBeta3DQGPerModel::addStatsOutputFiles(SharedSimulation spSim)
   {
   }

   void BoussinesqBeta3DQGPerModel::setInitialState(SharedSimulation spSim)
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
