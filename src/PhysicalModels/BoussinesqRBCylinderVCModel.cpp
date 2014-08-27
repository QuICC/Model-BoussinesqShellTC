/** 
 * @file BoussinesqRBCylinderVCModel.cpp
 * @brief Source of the Boussinesq Rayleigh-Benard cylinder (velocity-continuity formulation) model
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
#include "PhysicalModels/BoussinesqRBCylinderVCModel.hpp"

// Project includes
//
#include "Enums/FieldIds.hpp"
#include "IoVariable/StateFileReader.hpp"
#include "IoVariable/StateFileWriter.hpp"
#include "IoVariable/VisualizationFileWriter.hpp"
#include "IoTools/IdToHuman.hpp"
#include "Equations/Cylinder/Boussinesq/BoussinesqCylinderTransport.hpp"
#include "Equations/Cylinder/Boussinesq/BoussinesqCylinderVelocity.hpp"
#include "Generator/States/CylinderExactScalarState.hpp"
#include "Generator/States/CylinderExactVectorState.hpp"
#include "Generator/Visualizers/ScalarFieldVisualizer.hpp"
#include "Generator/Visualizers/VectorFieldVisualizer.hpp"
#include "PhysicalModels/PhysicalModelBase.hpp"

namespace GeoMHDiSCC {

   const std::string BoussinesqRBCylinderVCModel::PYMODULE = "boussinesq_rbcylinder_vc";

   const std::string BoussinesqRBCylinderVCModel::PYCLASS = "BoussinesqRBCylinderVC";

   void BoussinesqRBCylinderVCModel::addEquations(SharedSimulation spSim)
   {
      // Add transport equation
      spSim->addScalarEquation<Equations::BoussinesqCylinderTransport>();
      
      // Add Navier-Stokes equation
      spSim->addVectorEquation<Equations::BoussinesqCylinderVelocity>();
   }

   void BoussinesqRBCylinderVCModel::addStates(SharedStateGenerator spGen)
   {
      // Shared pointer to equation
      Equations::SharedCylinderExactScalarState spSExact;
      Equations::SharedCylinderExactVectorState spVExact;

      // Add temperature initial state generator
      spSExact = spGen->addScalarEquation<Equations::CylinderExactScalarState>();
      spSExact->setIdentity(PhysicalNames::TEMPERATURE);
      spSExact->setStateType(Equations::CylinderExactScalarState::CONSTANT);

      // Add temperature initial state generator
      spVExact = spGen->addVectorEquation<Equations::CylinderExactVectorState>();
      spVExact->setIdentity(PhysicalNames::VELOCITY);
      spVExact->setStateType(Equations::CylinderExactVectorState::CONSTANT);

      // Add output file
      IoVariable::SharedStateFileWriter spOut(new IoVariable::StateFileWriter(SchemeType::type(), SchemeType::isRegular()));
      spOut->expect(PhysicalNames::TEMPERATURE);
      spOut->expect(PhysicalNames::VELOCITY);
      spGen->addHdf5OutputFile(spOut);
   }

   void BoussinesqRBCylinderVCModel::addVisualizers(SharedVisualizationGenerator spVis)
   {
      // Shared pointer to basic field visualizer
      Equations::SharedScalarFieldVisualizer spField;

      // Add first field visualization
      spField = spVis->addScalarEquation<Equations::ScalarFieldVisualizer>();
      spField->setFields(true, false);
      spField->setIdentity(PhysicalNames::TEMPERATURE);

      // Add output file
      IoVariable::SharedVisualizationFileWriter spOut(new IoVariable::VisualizationFileWriter(SchemeType::type()));
      spOut->expect(PhysicalNames::TEMPERATURE);
      spVis->addHdf5OutputFile(spOut);
   }

   void BoussinesqRBCylinderVCModel::setVisualizationState(SharedVisualizationGenerator spVis)
   {
      // Create and add initial state file to IO
      IoVariable::SharedStateFileReader spIn(new IoVariable::StateFileReader("4Visu", SchemeType::type(), SchemeType::isRegular()));

      // Set expected fields
      spIn->expect(PhysicalNames::TEMPERATURE);

      // Set simulation state
      spVis->setInitialState(spIn);
   }

   void BoussinesqRBCylinderVCModel::addAsciiOutputFiles(SharedSimulation spSim)
   {
      // Add ASCII output file
      //pSim->addOutputFile(AN_ASCIIFILE);
      
      // Add ASCII output file
      //pSim->addOutputFile(AN_ASCIIFILE);
   }

   void BoussinesqRBCylinderVCModel::addHdf5OutputFiles(SharedSimulation spSim)
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

   void BoussinesqRBCylinderVCModel::setInitialState(SharedSimulation spSim)
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
