/** 
 * @file BoussinesqBeta3DQGModel.cpp
 * @brief Source of the Boussinesq Beta 3DQG model
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
#include "PhysicalModels/BoussinesqBeta3DQGModel.hpp"

// Project includes
//
#include "Enums/FieldIds.hpp"
#include "IoVariable/StateFileReader.hpp"
#include "IoVariable/StateFileWriter.hpp"
#include "IoVariable/VisualizationFileWriter.hpp"
#include "IoTools/IdToHuman.hpp"
#include "Equations/Asymptotics/Beta3DQG/Boussinesq/BoussinesqBeta3DQGStreamfunction.hpp"
#include "Equations/Asymptotics/Beta3DQG/Boussinesq/BoussinesqBeta3DQGVertical.hpp"
#include "Equations/Asymptotics/Beta3DQG/Boussinesq/BoussinesqBeta3DQGTransport.hpp"
#include "Equations/Asymptotics/Beta3DQG/Boussinesq/BoussinesqBeta3DQGVorticity.hpp"
#include "Generator/States/RandomScalarState.hpp"
#include "Generator/States/CartesianExactScalarState.hpp"
#include "Generator/Visualizers/ScalarFieldVisualizer.hpp"
#include "Generator/Visualizers/VorticityStreamVisualizer.hpp"
#include "PhysicalModels/PhysicalModelBase.hpp"

namespace GeoMHDiSCC {

   const std::string BoussinesqBeta3DQGModel::PYMODULE = "boussinesq_beta3dqg";

   const std::string BoussinesqBeta3DQGModel::PYCLASS = "BoussinesqBeta3DQG";

   void BoussinesqBeta3DQGModel::addEquations(SharedSimulation spSim)
   {
      // Add streamfunction equation
      spSim->addScalarEquation<Equations::BoussinesqBeta3DQGStreamfunction>();
      
      // Add vertical velocity equation
      spSim->addScalarEquation<Equations::BoussinesqBeta3DQGVertical>();
      
      // Add transport equation
      spSim->addScalarEquation<Equations::BoussinesqBeta3DQGTransport>();
      
      // Add vorticity computation
      spSim->addScalarEquation<Equations::BoussinesqBeta3DQGVorticity>();
   }

   void BoussinesqBeta3DQGModel::addStates(SharedStateGenerator spGen)
   {
      // Shared pointer to equation
      Equations::SharedRandomScalarState spRand;
      // Shared pointer to equation
      Equations::SharedCartesianExactScalarState spExact;

      // Add transport initial state generation equation
      spRand = spGen->addScalarEquation<Equations::RandomScalarState>();
      spRand->setIdentity(PhysicalNames::TEMPERATURE);
      spRand->setSpectrum(-0.1,0.1, 1e4, 1e4, 1e4);
      
      // Add streamfunction initial state generation equation
      spExact = spGen->addScalarEquation<Equations::CartesianExactScalarState>();
      spExact->setIdentity(PhysicalNames::STREAMFUNCTION);
      spExact->setStateType(Equations::CartesianExactScalarState::POLYSINPOLY);
      spExact->setModeOptions(1e0, 1.0, 1e0, 1.0, 1e0, 1.0);
      
      // Add vertical velocity initial state generation equation
      spExact = spGen->addScalarEquation<Equations::CartesianExactScalarState>();
      spExact->setIdentity(PhysicalNames::VELOCITYZ);
      spExact->setStateType(Equations::CartesianExactScalarState::POLYCOSPOLY);
      spExact->setModeOptions(1e0, 1.0, 1e0, 1.0, 1e0, 1.0);
      
//      // Add streamfunction initial state generation equation
//      spRand = spGen->addScalarEquation<Equations::RandomScalarState>();
//      spRand->setIdentity(PhysicalNames::STREAMFUNCTION);
//      spRand->setSpectrum(-0.1,0.1, 1e4, 1e4, 1e4);
//      
//      // Add vertical velocity initial state generation equation
//      spRand = spGen->addScalarEquation<Equations::RandomScalarState>();
//      spRand->setIdentity(PhysicalNames::VELOCITYZ);
//      spRand->setSpectrum(-0.1,0.1, 1e4, 1e4, 1e4);

      // Add output file
      IoVariable::SharedStateFileWriter spOut(new IoVariable::StateFileWriter(SchemeType::type(), SchemeType::isRegular()));
      spOut->expect(PhysicalNames::TEMPERATURE);
      spOut->expect(PhysicalNames::STREAMFUNCTION);
      spOut->expect(PhysicalNames::VELOCITYZ);
      spGen->addHdf5OutputFile(spOut);
   }

   void BoussinesqBeta3DQGModel::addVisualizers(SharedVisualizationGenerator spVis)
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
      
      // Add vorticity field visualization
      Equations::SharedVorticityStreamVisualizer spVort;
      spVort = spVis->addScalarEquation<Equations::VorticityStreamVisualizer>();
      spVort->setFields(true, true);

      // Add output file
      IoVariable::SharedVisualizationFileWriter spOut(new IoVariable::VisualizationFileWriter(SchemeType::type()));
      spOut->expect(PhysicalNames::TEMPERATURE);
      spOut->expect(PhysicalNames::STREAMFUNCTION);
      spOut->expect(PhysicalNames::VELOCITYZ);
      spOut->expect(PhysicalNames::VORTICITY);
      spVis->addHdf5OutputFile(spOut);
   }

   void BoussinesqBeta3DQGModel::setVisualizationState(SharedVisualizationGenerator spVis)
   {
      // Create and add initial state file to IO
      IoVariable::SharedStateFileReader spIn(new IoVariable::StateFileReader("4Visu", SchemeType::type(), SchemeType::isRegular()));

      // Set expected fields
      spIn->expect(PhysicalNames::TEMPERATURE);
      spIn->expect(PhysicalNames::STREAMFUNCTION);
      spIn->expect(PhysicalNames::VELOCITYZ);

      // Set simulation state
      spVis->setInitialState(spIn);
   }

   void BoussinesqBeta3DQGModel::addAsciiOutputFiles(SharedSimulation spSim)
   {
      // Add ASCII output file
      //pSim->addOutputFile(AN_ASCIIFILE);
      
      // Add ASCII output file
      //pSim->addOutputFile(AN_ASCIIFILE);
   }

   void BoussinesqBeta3DQGModel::addHdf5OutputFiles(SharedSimulation spSim)
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

   void BoussinesqBeta3DQGModel::setInitialState(SharedSimulation spSim)
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
