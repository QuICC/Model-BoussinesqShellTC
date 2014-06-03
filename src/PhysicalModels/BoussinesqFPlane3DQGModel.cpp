/** 
 * @file BoussinesqFPlane3DQGModel.cpp
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
#include "PhysicalModels/BoussinesqFPlane3DQGModel.hpp"

// Project includes
//
#include "Enums/FieldIds.hpp"
#include "IoVariable/StateFileReader.hpp"
#include "IoVariable/StateFileWriter.hpp"
#include "IoVariable/VisualizationFileWriter.hpp"
#include "IoTools/IdToHuman.hpp"
#include "Equations/Asymptotics/FPlane3DQG/Boussinesq/BoussinesqFPlane3DQGMeanHeat.hpp"
#include "Equations/Asymptotics/FPlane3DQG/Boussinesq/BoussinesqFPlane3DQGStreamfunction.hpp"
#include "Equations/Asymptotics/FPlane3DQG/Boussinesq/BoussinesqFPlane3DQGTransport.hpp"
#include "Equations/Asymptotics/FPlane3DQG/Boussinesq/BoussinesqFPlane3DQGVertical.hpp"
#include "Equations/Asymptotics/FPlane3DQG/Boussinesq/BoussinesqFPlane3DQGVorticity.hpp"
#include "Generator/States/RandomScalarState.hpp"
#include "Generator/States/CartesianExactScalarState.hpp"
#include "Generator/Visualizers/ScalarFieldVisualizer.hpp"
#include "PhysicalModels/PhysicalModelBase.hpp"

namespace GeoMHDiSCC {

   const std::string BoussinesqFPlane3DQGModel::PYNAME = "boussinesq_fplane3dqg_model";

   void BoussinesqFPlane3DQGModel::addEquations(SharedSimulation spSim)
   {
      // Add streamfunction equation
      spSim->addScalarEquation<Equations::BoussinesqFPlane3DQGStreamfunction>(BoussinesqFPlane3DQGModel::PYNAME);
      
      // Add vertical velocity equation
      spSim->addScalarEquation<Equations::BoussinesqFPlane3DQGVertical>(BoussinesqFPlane3DQGModel::PYNAME);
      
      // Add transport equation
      spSim->addScalarEquation<Equations::BoussinesqFPlane3DQGTransport>(BoussinesqFPlane3DQGModel::PYNAME);
      
      // Add vorticity computation
      spSim->addScalarEquation<Equations::BoussinesqFPlane3DQGVorticity>(BoussinesqFPlane3DQGModel::PYNAME);
      
      // Add mean heat computation
      spSim->addScalarEquation<Equations::BoussinesqFPlane3DQGMeanHeat>(BoussinesqFPlane3DQGModel::PYNAME);
   }

   void BoussinesqFPlane3DQGModel::addStates(SharedStateGenerator spGen)
   {
      // Shared pointer to equation
      Equations::SharedCartesianExactScalarState spExact;

      // Add transport initial state generation equation
      spExact = spGen->addScalarEquation<Equations::CartesianExactScalarState>(BoussinesqFPlane3DQGModel::PYNAME);
      spExact->setIdentity(PhysicalNames::TEMPERATURE);
      spExact->setStateType(Equations::CartesianExactScalarState::POLYCOSCOS);
      spExact->setModeOptions(1e0, 1.0, 1e0, 1.0, 1e0, 1.0);
      
      // Add streamfunction initial state generation equation
      spExact = spGen->addScalarEquation<Equations::CartesianExactScalarState>(BoussinesqFPlane3DQGModel::PYNAME);
      spExact->setIdentity(PhysicalNames::STREAMFUNCTION);
      spExact->setStateType(Equations::CartesianExactScalarState::POLYSINCOS);
      spExact->setModeOptions(1e0, 2.0, 1e0, 1.0, 1e0, 1.0);
      
      // Add vertical velocity initial state generation equation
      spExact = spGen->addScalarEquation<Equations::CartesianExactScalarState>(BoussinesqFPlane3DQGModel::PYNAME);
      spExact->setIdentity(PhysicalNames::VELOCITYZ);
      spExact->setStateType(Equations::CartesianExactScalarState::POLYCOSSIN);
      spExact->setModeOptions(1e0, 1.0, 1e0, 1.0, 1e0, 2.0);

      // Add output file
      IoVariable::SharedStateFileWriter spOut(new IoVariable::StateFileWriter(SchemeType::type(), SchemeType::isRegular()));
      spOut->expect(PhysicalNames::TEMPERATURE);
      spOut->expect(PhysicalNames::STREAMFUNCTION);
      spOut->expect(PhysicalNames::VELOCITYZ);
      spGen->addOutputFile(spOut);
   }

   void BoussinesqFPlane3DQGModel::addVisualizers(SharedVisualizationGenerator spVis)
   {
      // Shared pointer to basic field visualizer
      Equations::SharedScalarFieldVisualizer spField;

      // Add transport field visualization
      spField = spVis->addScalarEquation<Equations::ScalarFieldVisualizer>(BoussinesqFPlane3DQGModel::PYNAME);
      spField->setFields(true, true);
      spField->setIdentity(PhysicalNames::TEMPERATURE);
      
      // Add streamfunction field visualization
      spField = spVis->addScalarEquation<Equations::ScalarFieldVisualizer>(BoussinesqFPlane3DQGModel::PYNAME);
      spField->setFields(true, true);
      spField->setIdentity(PhysicalNames::STREAMFUNCTION);
      
      // Add vertical velocity field visualization
      spField = spVis->addScalarEquation<Equations::ScalarFieldVisualizer>(BoussinesqFPlane3DQGModel::PYNAME);
      spField->setFields(true, true);
      spField->setIdentity(PhysicalNames::VELOCITYZ);
      
      // Add vertical velocity field visualization
      spField = spVis->addScalarEquation<Equations::ScalarFieldVisualizer>(BoussinesqFPlane3DQGModel::PYNAME);
      spField->setFields(true, true);
      spField->setIdentity(PhysicalNames::MEANTEMPERATURE);

      // Add output file
      IoVariable::SharedVisualizationFileWriter spOut(new IoVariable::VisualizationFileWriter(SchemeType::type()));
      spOut->expect(PhysicalNames::TEMPERATURE);
      spOut->expect(PhysicalNames::STREAMFUNCTION);
      spOut->expect(PhysicalNames::VELOCITYZ);
      spOut->expect(PhysicalNames::MEANTEMPERATURE);
      spVis->addOutputFile(spOut);
   }

   void BoussinesqFPlane3DQGModel::setVisualizationState(SharedVisualizationGenerator spVis)
   {
      // Create and add initial state file to IO
      IoVariable::SharedStateFileReader spIn(new IoVariable::StateFileReader("4Visu", SchemeType::type(), SchemeType::isRegular()));

      // Set expected fields
      spIn->expect(PhysicalNames::TEMPERATURE);
      spIn->expect(PhysicalNames::STREAMFUNCTION);
      spIn->expect(PhysicalNames::VELOCITYZ);
      spIn->expect(PhysicalNames::MEANTEMPERATURE);

      // Set simulation state
      spVis->setInitialState(spIn);
   }

   void BoussinesqFPlane3DQGModel::addAsciiOutputFiles(SharedSimulation spSim)
   {
      // Add ASCII output file
      //pSim->addOutputFile(AN_ASCIIFILE);
      
      // Add ASCII output file
      //pSim->addOutputFile(AN_ASCIIFILE);
   }

   void BoussinesqFPlane3DQGModel::addHdf5OutputFiles(SharedSimulation spSim)
   {
      // Field IDs iterator
      std::vector<GeoMHDiSCC::PhysicalNames::Id>::const_iterator  it;
      std::vector<GeoMHDiSCC::PhysicalNames::Id> ids = PhysicalModelBase::fieldIds(BoussinesqFPlane3DQGModel::PYNAME);

      // Create and add state file to IO
      IoVariable::SharedStateFileWriter spState(new IoVariable::StateFileWriter(SchemeType::type(), SchemeType::isRegular()));
      for(it = ids.begin(); it != ids.end(); ++it)
      {
         spState->expect(*it);
      }

      // Add mean temperature to ouput file
      spState->expect(PhysicalNames::MEANTEMPERATURE);

      spSim->addOutputFile(spState);
   }

   void BoussinesqFPlane3DQGModel::setInitialState(SharedSimulation spSim)
   {
      // Field IDs iterator
      std::vector<GeoMHDiSCC::PhysicalNames::Id>::const_iterator  it;
      std::vector<GeoMHDiSCC::PhysicalNames::Id> ids = PhysicalModelBase::fieldIds(BoussinesqFPlane3DQGModel::PYNAME);

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
