/** 
 * @file BoussinesqShellModel.cpp
 * @brief Source of the Boussinesq shell physical model
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
#include "PhysicalModels/BoussinesqShellModel.hpp"

// Project includes
//
#include "Enums/FieldIds.hpp"
#include "IoVariable/StateFileReader.hpp"
#include "IoVariable/StateFileWriter.hpp"
#include "IoVariable/VisualizationFileWriter.hpp"
#include "IoTools/IdToHuman.hpp"
#include "Equations/Shell/Boussinesq/BoussinesqShellTransport.hpp"
#include "Equations/Shell/Boussinesq/BoussinesqShellVelocity.hpp"
#include "Generator/States/ShellExactScalarState.hpp"
#include "Generator/States/ShellExactVectorState.hpp"
#include "Generator/Visualizers/ScalarFieldVisualizer.hpp"
#include "Generator/Visualizers/VectorFieldVisualizer.hpp"
#include "PhysicalModels/PhysicalModelBase.hpp"

namespace GeoMHDiSCC {

   const std::string BoussinesqShellModel::PYNAME = "boussinesq_shell_model";

   void BoussinesqShellModel::addEquations(SharedSimulation spSim)
   {
      // Add transport equation
      spSim->addScalarEquation<Equations::BoussinesqShellTransport>(BoussinesqShellModel::PYNAME);
      
      // Add Navier-Stokes equation
      spSim->addVectorEquation<Equations::BoussinesqShellVelocity>(BoussinesqShellModel::PYNAME);
   }

   void BoussinesqShellModel::addStates(SharedStateGenerator spGen)
   {
      // Shared pointer to equation
      Equations::SharedShellExactScalarState spSExact;
      Equations::SharedShellExactVectorState spVExact;

      // Add temperature initial state generator
      spSExact = spGen->addScalarEquation<Equations::ShellExactScalarState>();
      spSExact->setIdentity(PhysicalNames::TEMPERATURE);
      spSExact->setStateType(Equations::ShellExactScalarState::HARMONIC);
      std::vector<std::tr1::tuple<int,int,MHDComplex> > tSH;
      tSH.push_back(std::tr1::make_tuple(23,9,MHDComplex(1,1)));
      spSExact->setHarmonicOptions(tSH);

      // Add temperature initial state generator
      spVExact = spGen->addVectorEquation<Equations::ShellExactVectorState>();
      spVExact->setIdentity(PhysicalNames::VELOCITY);
      spVExact->setStateType(Equations::ShellExactVectorState::CONSTANT);

      // Add output file
      IoVariable::SharedStateFileWriter spOut(new IoVariable::StateFileWriter(SchemeType::type(), SchemeType::isRegular()));
      spOut->expect(PhysicalNames::TEMPERATURE);
      spOut->expect(PhysicalNames::VELOCITY);
      spGen->addOutputFile(spOut);
   }

   void BoussinesqShellModel::addVisualizers(SharedVisualizationGenerator spVis)
   {
      // Shared pointer to basic field visualizer
      Equations::SharedScalarFieldVisualizer spSField;
      Equations::SharedVectorFieldVisualizer spVField;

      // Add temperature field visualization
      spSField = spVis->addScalarEquation<Equations::ScalarFieldVisualizer>();
      spSField->setFields(true, false);
      spSField->setIdentity(PhysicalNames::TEMPERATURE);

      // Add velocity field visualization
      spVField = spVis->addVectorEquation<Equations::VectorFieldVisualizer>();
      spVField->setFields(true, false);
      spVField->setIdentity(PhysicalNames::VELOCITY);

      // Add output file
      IoVariable::SharedVisualizationFileWriter spOut(new IoVariable::VisualizationFileWriter(SchemeType::type()));
      spOut->expect(PhysicalNames::TEMPERATURE);
      spOut->expect(PhysicalNames::VELOCITY);
      spVis->addOutputFile(spOut);
   }

   void BoussinesqShellModel::setVisualizationState(SharedVisualizationGenerator spVis)
   {
      // Create and add initial state file to IO
      IoVariable::SharedStateFileReader spIn(new IoVariable::StateFileReader("4Visu", SchemeType::type(), SchemeType::isRegular()));

      // Set expected fields
      spIn->expect(PhysicalNames::TEMPERATURE);

      // Set simulation state
      spVis->setInitialState(spIn);
   }

   void BoussinesqShellModel::addAsciiOutputFiles(SharedSimulation spSim)
   {
      // Add ASCII output file
      //pSim->addOutputFile(AN_ASCIIFILE);
      
      // Add ASCII output file
      //pSim->addOutputFile(AN_ASCIIFILE);
   }

   void BoussinesqShellModel::addHdf5OutputFiles(SharedSimulation spSim)
   {
      // Field IDs iterator
      std::vector<GeoMHDiSCC::PhysicalNames::Id>::const_iterator  it;
      std::vector<GeoMHDiSCC::PhysicalNames::Id> ids = PhysicalModelBase::fieldIds(BoussinesqShellModel::PYNAME);

      // Create and add state file to IO
      IoVariable::SharedStateFileWriter spState(new IoVariable::StateFileWriter(SchemeType::type(), SchemeType::isRegular()));
      for(it = ids.begin(); it != ids.end(); ++it)
      {
         spState->expect(*it);
      }
      spSim->addOutputFile(spState);
   }

   void BoussinesqShellModel::setInitialState(SharedSimulation spSim)
   {
      // Field IDs iterator
      std::vector<GeoMHDiSCC::PhysicalNames::Id>::const_iterator  it;
      std::vector<GeoMHDiSCC::PhysicalNames::Id> ids = PhysicalModelBase::fieldIds(BoussinesqShellModel::PYNAME);

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
