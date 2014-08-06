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
#include "IoVariable/NusseltWriter.hpp"
#include "IoVariable/VisualizationFileWriter.hpp"
#include "IoTools/IdToHuman.hpp"
#include "Equations/Asymptotics/Beta3DQG/Boussinesq/BoussinesqBeta3DQGPerStreamfunction.hpp"
#include "Equations/Asymptotics/Beta3DQG/Boussinesq/BoussinesqBeta3DQGPerVelocityZ.hpp"
#include "Equations/Asymptotics/Beta3DQG/Boussinesq/BoussinesqBeta3DQGPerTransport.hpp"
#include "Equations/Asymptotics/Beta3DQG/Boussinesq/BoussinesqBeta3DQGPerVorticityZ.hpp"
#include "Generator/States/RandomScalarState.hpp"
#include "Generator/States/CartesianExactScalarState.hpp"
#include "Generator/Visualizers/ScalarFieldVisualizer.hpp"
#include "PhysicalModels/PhysicalModelBase.hpp"

namespace GeoMHDiSCC {

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
         spExact->setStateType(Equations::CartesianExactScalarState::POLYCOSCOS);
         spExact->setModeOptions(1e0, 0.0, 1e0, 0.0, 1e0, 0.0);

         // Add streamfunction initial state generation equation
         spExact = spGen->addScalarEquation<Equations::CartesianExactScalarState>();
         spExact->setIdentity(PhysicalNames::STREAMFUNCTION);
         spExact->setStateType(Equations::CartesianExactScalarState::POLYCOSCOS);
         spExact->setModeOptions(1e0, 1.0, 1e0, 0.0, 1e0, 0.0);

         // Add vertical velocity initial state generation equation
         spExact = spGen->addScalarEquation<Equations::CartesianExactScalarState>();
         spExact->setIdentity(PhysicalNames::VELOCITYZ);
         spExact->setStateType(Equations::CartesianExactScalarState::POLYCOSCOS);
         spExact->setModeOptions(1e0, 2.0, 1e0, 0.0, 1e0, 0.0);

         // Add vertical vorticity initial state generation equation
         spExact = spGen->addScalarEquation<Equations::CartesianExactScalarState>();
         spExact->setIdentity(PhysicalNames::VORTICITYZ);
         spExact->setStateType(Equations::CartesianExactScalarState::POLYCOSCOS);
         spExact->setModeOptions(1e0, 2.0, 1e0, 0.0, 1e0, 0.0);

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
      spField->setFields(true, true);
      spField->setIdentity(PhysicalNames::TEMPERATURE);
      
      // Add streamfunction field visualization
      spField = spVis->addScalarEquation<Equations::ScalarFieldVisualizer>();
      spField->setFields(true, true);
      spField->setIdentity(PhysicalNames::STREAMFUNCTION);
      
      // Add vertical velocity field visualization
      spField = spVis->addScalarEquation<Equations::ScalarFieldVisualizer>();
      spField->setFields(true, true);
      spField->setIdentity(PhysicalNames::VELOCITYZ);
      
      // Add vertical velocity field visualization
      spField = spVis->addScalarEquation<Equations::ScalarFieldVisualizer>();
      spField->setFields(true, true);
      spField->setIdentity(PhysicalNames::VORTICITYZ);

      // Add output file
      IoVariable::SharedVisualizationFileWriter spOut(new IoVariable::VisualizationFileWriter(SchemeType::type()));
      spOut->expect(PhysicalNames::TEMPERATURE);
      spOut->expect(PhysicalNames::STREAMFUNCTION);
      spOut->expect(PhysicalNames::VELOCITYZ);
      spOut->expect(PhysicalNames::VORTICITYZ);
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

      // Set simulation state
      spVis->setInitialState(spIn);
   }

   void BoussinesqBeta3DQGPerModel::addAsciiOutputFiles(SharedSimulation spSim)
   {
      // Add ASCII output file
      //pSim->addOutputFile(AN_ASCIIFILE);
   }

   void BoussinesqBeta3DQGPerModel::addHdf5OutputFiles(SharedSimulation spSim)
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

   void BoussinesqBeta3DQGPerModel::setInitialState(SharedSimulation spSim)
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
