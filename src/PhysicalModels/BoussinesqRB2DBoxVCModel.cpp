/** 
 * @file BoussinesqRB2DBoxVCModel.cpp
 * @brief Source of the Boussinesq Rayleigh-Benard 2D box (velocity-continuity formulation) model
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
#include "PhysicalModels/BoussinesqRB2DBoxVCModel.hpp"

// Project includes
//
#include "Enums/FieldIds.hpp"
#include "IoVariable/StateFileReader.hpp"
#include "IoVariable/StateFileWriter.hpp"
#include "IoVariable/VisualizationFileWriter.hpp"
#include "IoTools/IdToHuman.hpp"
#include "Equations/Box/Boussinesq/BoussinesqBoxTransport.hpp"
#include "Equations/Box/Boussinesq/BoussinesqBoxVelocity.hpp"
#include "Generator/States/CartesianExactScalarState.hpp"
#include "Generator/States/CartesianExactVectorState.hpp"
#include "Generator/Visualizers/ScalarFieldVisualizer.hpp"
#include "Generator/Visualizers/VectorFieldVisualizer.hpp"
#include "PhysicalModels/PhysicalModelBase.hpp"

namespace GeoMHDiSCC {

   const std::string BoussinesqRB2DBoxVCModel::PYMODULE = "boussinesq_rb2dbox_vc";

   const std::string BoussinesqRB2DBoxVCModel::PYCLASS = "BoussinesqRB2DBoxVC";

   void BoussinesqRB2DBoxVCModel::addEquations(SharedSimulation spSim)
   {
      // Add transport equation
      spSim->addScalarEquation<Equations::BoussinesqBoxTransport>();
      
      // Add Navier-Stokes equation
      spSim->addScalarEquation<Equations::BoussinesqBoxVelocityX>();
      spSim->addScalarEquation<Equations::BoussinesqBoxVelocityY>();
      spSim->addScalarEquation<Equations::BoussinesqBoxVelocityZ>();
   }

   void BoussinesqRB2DBoxVCModel::addStates(SharedStateGenerator spGen)
   {
      // Shared pointer to equation
      Equations::SharedCartesianExactScalarState spSExact;

      // Add temperature initial state generator
      spSExact = spGen->addScalarEquation<Equations::CartesianExactScalarState>();
      spSExact->setIdentity(PhysicalNames::TEMPERATURE);
      spSExact->setStateType(Equations::CartesianExactScalarState::CONSTANT);

      // Add velocity initial state generators
      spSExact = spGen->addScalarEquation<Equations::CartesianExactScalarState>();
      spSExact->setIdentity(PhysicalNames::VELOCITYX);
      spSExact->setStateType(Equations::CartesianExactScalarState::CONSTANT);
      spSExact = spGen->addScalarEquation<Equations::CartesianExactScalarState>();
      spSExact->setIdentity(PhysicalNames::VELOCITYY);
      spSExact->setStateType(Equations::CartesianExactScalarState::CONSTANT);
      spSExact = spGen->addScalarEquation<Equations::CartesianExactScalarState>();
      spSExact->setIdentity(PhysicalNames::VELOCITYZ);
      spSExact->setStateType(Equations::CartesianExactScalarState::CONSTANT);

      // Add output file
      IoVariable::SharedStateFileWriter spOut(new IoVariable::StateFileWriter(SchemeType::type(), SchemeType::isRegular()));
      spOut->expect(PhysicalNames::TEMPERATURE);
      spOut->expect(PhysicalNames::VELOCITYX);
      spOut->expect(PhysicalNames::VELOCITYY);
      spOut->expect(PhysicalNames::VELOCITYZ);
      spGen->addHdf5OutputFile(spOut);
   }

   void BoussinesqRB2DBoxVCModel::addVisualizers(SharedVisualizationGenerator spVis)
   {
      // Shared pointer to basic field visualizer
      Equations::SharedScalarFieldVisualizer spField;

      // Add temperature field visualization
      spField = spVis->addScalarEquation<Equations::ScalarFieldVisualizer>();
      spField->setFields(true, false);
      spField->setIdentity(PhysicalNames::TEMPERATURE);

      // Add velocity fields visualization
      spField = spVis->addScalarEquation<Equations::ScalarFieldVisualizer>();
      spField->setFields(true, false);
      spField->setIdentity(PhysicalNames::VELOCITYX);
      spField = spVis->addScalarEquation<Equations::ScalarFieldVisualizer>();
      spField->setFields(true, false);
      spField->setIdentity(PhysicalNames::VELOCITYY);
      spField = spVis->addScalarEquation<Equations::ScalarFieldVisualizer>();
      spField->setFields(true, false);
      spField->setIdentity(PhysicalNames::VELOCITYZ);

      // Add output file
      IoVariable::SharedVisualizationFileWriter spOut(new IoVariable::VisualizationFileWriter(SchemeType::type()));
      spOut->expect(PhysicalNames::TEMPERATURE);
      spVis->addHdf5OutputFile(spOut);
   }

   void BoussinesqRB2DBoxVCModel::setVisualizationState(SharedVisualizationGenerator spVis)
   {
      // Create and add initial state file to IO
      IoVariable::SharedStateFileReader spIn(new IoVariable::StateFileReader("4Visu", SchemeType::type(), SchemeType::isRegular()));

      // Set expected fields
      spIn->expect(PhysicalNames::TEMPERATURE);
      spIn->expect(PhysicalNames::VELOCITYX);
      spIn->expect(PhysicalNames::VELOCITYY);
      spIn->expect(PhysicalNames::VELOCITYZ);

      // Set simulation state
      spVis->setInitialState(spIn);
   }

   void BoussinesqRB2DBoxVCModel::addAsciiOutputFiles(SharedSimulation spSim)
   {
      // Add ASCII output file
      //pSim->addOutputFile(AN_ASCIIFILE);
      
      // Add ASCII output file
      //pSim->addOutputFile(AN_ASCIIFILE);
   }

   void BoussinesqRB2DBoxVCModel::addHdf5OutputFiles(SharedSimulation spSim)
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

   void BoussinesqRB2DBoxVCModel::setInitialState(SharedSimulation spSim)
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
