/** 
 * @file BoussinesqBetaCylGModel.cpp
 * @brief Source of the Boussinesq beta 3DQG model with cylindrical gravity
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
#include "PhysicalModels/BoussinesqBetaCylGModel.hpp"

// Project includes
//
#include "Enums/FieldIds.hpp"
#include "BoundaryCondition/BoundaryCondition.hpp"
#include "IoVariable/StateFileReader.hpp"
#include "IoVariable/StateFileWriter.hpp"
#include "IoVariable/VisualizationFileWriter.hpp"
#include "IoTools/IdToHuman.hpp"
#include "Equations/Asymptotics/Beta3DQG/Boussinesq/BoussinesqBetaCylGStreamfunction.hpp"
#include "Equations/Asymptotics/Beta3DQG/Boussinesq/BoussinesqBetaCylGVertical.hpp"
#include "Equations/Asymptotics/Beta3DQG/Boussinesq/BoussinesqBetaCylGTransport.hpp"
#include "Equations/Asymptotics/Beta3DQG/Boussinesq/BoussinesqBetaCylGVorticity.hpp"
#include "Generator/States/RandomScalarState.hpp"
#include "Generator/States/ExactScalarState.hpp"
#include "Generator/Visualizers/ScalarFieldVisualizer.hpp"
#include "Generator/Visualizers/VorticityStreamVisualizer.hpp"
#include "PhysicalModels/PhysicalModelBase.hpp"

namespace GeoMHDiSCC {

   const std::string BoussinesqBetaCylGModel::PYNAME = "boussinesq_beta_cylg_model";

   void BoussinesqBetaCylGModel::addEquations(SharedSimulation spSim)
   {
      // Add transport equation
      spSim->addScalarEquation<Equations::BoussinesqBetaCylGTransport>(BoussinesqBetaCylGModel::PYNAME);
      
      // Add streamfunction equation
      spSim->addScalarEquation<Equations::BoussinesqBetaCylGStreamfunction>(BoussinesqBetaCylGModel::PYNAME);
      
      // Add vertical velocity equation
      spSim->addScalarEquation<Equations::BoussinesqBetaCylGVertical>(BoussinesqBetaCylGModel::PYNAME);
      
      // Add vorticity computation
      spSim->addScalarEquation<Equations::BoussinesqBetaCylGVorticity>(BoussinesqBetaCylGModel::PYNAME);
   }

   void BoussinesqBetaCylGModel::addStates(SharedStateGenerator spGen)
   {
      // Shared pointer to equation
      Equations::SharedRandomScalarState spRand;
      // Shared pointer to equation
      Equations::SharedExactScalarState spExact;

      // Add transport initial state generation equation
      spRand = spGen->addScalarEquation<Equations::RandomScalarState>();
      spRand->setIdentity(PhysicalNames::TEMPERATURE);
      spRand->setSpectrum(-0.1,0.1, 1e4, 1e4, 1e4);
      
      // Add streamfunction initial state generation equation
      spExact = spGen->addScalarEquation<Equations::ExactScalarState>();
      spExact->setIdentity(PhysicalNames::STREAMFUNCTION);
      spExact->setStateType(Equations::ExactScalarState::SINECOSINE);
      spExact->setSineOptions(1.0, 3, 1.0, 7.0);
      
      // Add vertical velocity initial state generation equation
      spExact = spGen->addScalarEquation<Equations::ExactScalarState>();
      spExact->setIdentity(PhysicalNames::VELOCITYZ);
      spExact->setStateType(Equations::ExactScalarState::SINESINE);
      spExact->setSineOptions(1.0, 3, 1.0, 6.0);
      
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
      spGen->addOutputFile(spOut);
   }

   void BoussinesqBetaCylGModel::addVisualizers(SharedVisualizationGenerator spVis)
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
      spVis->addOutputFile(spOut);
   }

   void BoussinesqBetaCylGModel::setVisualizationState(SharedVisualizationGenerator spVis)
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

   void BoussinesqBetaCylGModel::addAsciiOutputFiles(SharedSimulation spSim)
   {
      // Add ASCII output file
      //pSim->addOutputFile(AN_ASCIIFILE);
      
      // Add ASCII output file
      //pSim->addOutputFile(AN_ASCIIFILE);
   }

   void BoussinesqBetaCylGModel::addHdf5OutputFiles(SharedSimulation spSim)
   {
      // Field IDs iterator
      std::vector<GeoMHDiSCC::PhysicalNames::Id>::const_iterator  it;
      std::vector<GeoMHDiSCC::PhysicalNames::Id> ids = PhysicalModelBase::fieldIds(BoussinesqBetaCylGModel::PYNAME);

      // Create and add state file to IO
      IoVariable::SharedStateFileWriter spState(new IoVariable::StateFileWriter(SchemeType::type(), SchemeType::isRegular()));
      for(it = ids.begin(); it != ids.end(); ++it)
      {
         spState->expect(*it);
      }
      spSim->addOutputFile(spState);
   }

   void BoussinesqBetaCylGModel::setInitialState(SharedSimulation spSim)
   {
      // Field IDs iterator
      std::vector<GeoMHDiSCC::PhysicalNames::Id>::const_iterator  it;
      std::vector<GeoMHDiSCC::PhysicalNames::Id> ids = PhysicalModelBase::fieldIds(BoussinesqBetaCylGModel::PYNAME);

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
