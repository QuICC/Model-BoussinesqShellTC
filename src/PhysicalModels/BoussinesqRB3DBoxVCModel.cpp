/** 
 * @file BoussinesqRB3DBoxVCModel.cpp
 * @brief Source of the Boussinesq Rayleigh-Benard 3D box (velocity-continuity formulation) model
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
#include "PhysicalModels/BoussinesqRB3DBoxVCModel.hpp"

// Project includes
//
#include "Enums/FieldIds.hpp"
#include "IoVariable/StateFileReader.hpp"
#include "IoVariable/StateFileWriter.hpp"
#include "IoVariable/ContinuityWriter.hpp"
#include "IoVariable/NusseltCubicWriter.hpp"
#include "IoVariable/VisualizationFileWriter.hpp"
#include "IoTools/IdToHuman.hpp"
#include "Equations/Box/Boussinesq/BoussinesqRB3DBoxVCTransport.hpp"
#include "Equations/Box/Boussinesq/BoussinesqRB3DBoxVCMomentum.hpp"
#include "Equations/Box/Boussinesq/BoussinesqRB3DBoxVCContinuity.hpp"
#include "Generator/States/RandomScalarState.hpp"
#include "Generator/States/CartesianExactScalarState.hpp"
#include "Generator/Visualizers/ScalarFieldVisualizer.hpp"
#include "PhysicalModels/PhysicalModelBase.hpp"

namespace GeoMHDiSCC {

   const std::string BoussinesqRB3DBoxVCModel::PYMODULE = "boussinesq_rb3dbox_vc";

   const std::string BoussinesqRB3DBoxVCModel::PYCLASS = "BoussinesqRB3DBoxVC";

   void BoussinesqRB3DBoxVCModel::addEquations(SharedSimulation spSim)
   {
      // Add transport equation
      spSim->addScalarEquation<Equations::BoussinesqRB3DBoxVCTransport>();
      
      // Add Navier-Stokes equation (X,Y,Z components)
      spSim->addVectorEquation<Equations::BoussinesqRB3DBoxVCMomentum>();

      // Add continuity equation
      spSim->addScalarEquation<Equations::BoussinesqRB3DBoxVCContinuity>();
   }

   void BoussinesqRB3DBoxVCModel::addStates(SharedStateGenerator spGen)
   {
      // Generate "exact" solutions (trigonometric or monomial)
      if(false)
      {
         // Shared pointer to equation
         Equations::SharedCartesianExactScalarState spExact;

         // Add scalar exact initial state generator
         spExact = spGen->addScalarEquation<Equations::CartesianExactScalarState>();
         spExact->setIdentity(PhysicalNames::VELOCITYX);
         spExact->setStateType(Equations::CartesianExactScalarState::POLYPOLYPOLY);
         spExact->setModeOptions(1e0, 1.0, 1e0, 0.0, 1e0, 0.0);

         // Add scalar exact initial state generator
         spExact = spGen->addScalarEquation<Equations::CartesianExactScalarState>();
         spExact->setIdentity(PhysicalNames::VELOCITYY);
         spExact->setStateType(Equations::CartesianExactScalarState::POLYPOLYPOLY);
         spExact->setModeOptions(1e0, 0.0, 1e0, 1.0, 1e0, 0.0);

         // Add scalar exact initial state generator
         spExact = spGen->addScalarEquation<Equations::CartesianExactScalarState>();
         spExact->setIdentity(PhysicalNames::VELOCITYZ);
         spExact->setStateType(Equations::CartesianExactScalarState::POLYPOLYPOLY);
         spExact->setModeOptions(1e0, 0.0, 1e0, 0.0, 1e0, 1.0);

         // Add scalar exact initial state generator
         spExact = spGen->addScalarEquation<Equations::CartesianExactScalarState>();
         spExact->setIdentity(PhysicalNames::TEMPERATURE);
         spExact->setStateType(Equations::CartesianExactScalarState::POLYPOLYPOLY);
         spExact->setModeOptions(1e1, 1.0, 1e1, 2.0, 1e1, 3.0);

      // Generate random spectrum
      } else
      {
         // Shared pointer to random initial state equation
         Equations::SharedRandomScalarState spRand;

         // Add scalar random initial state generator 
         spRand = spGen->addScalarEquation<Equations::RandomScalarState>();
         spRand->setIdentity(PhysicalNames::VELOCITYX);
         spRand->setSpectrum(-1e-4, 1e-4, 1e4, 1e4, 1e4);

         // Add scalar random initial state generator 
         spRand = spGen->addScalarEquation<Equations::RandomScalarState>();
         spRand->setIdentity(PhysicalNames::VELOCITYY);
         spRand->setSpectrum(-1e-4, 1e-4, 1e4, 1e4, 1e4);

         // Add scalar random initial state generator 
         spRand = spGen->addScalarEquation<Equations::RandomScalarState>();
         spRand->setIdentity(PhysicalNames::VELOCITYZ);
         spRand->setSpectrum(-1e-4, 1e-4, 1e4, 1e4, 1e4);

         // Add scalar random initial state generator
         spRand = spGen->addScalarEquation<Equations::RandomScalarState>();
         spRand->setIdentity(PhysicalNames::TEMPERATURE);
         spRand->setSpectrum(-1e-4, 1e-4, 1e4, 1e4, 1e4);
      }

      // Add output file
      IoVariable::SharedStateFileWriter spOut(new IoVariable::StateFileWriter(SchemeType::type(), SchemeType::isRegular()));
      spOut->expect(PhysicalNames::VELOCITY);
      spOut->expect(PhysicalNames::TEMPERATURE);
      spGen->addHdf5OutputFile(spOut);
   }

   void BoussinesqRB3DBoxVCModel::addVisualizers(SharedVisualizationGenerator spVis)
   {
      // Shared pointer to basic field visualizer
      Equations::SharedScalarFieldVisualizer spField;

      // Add temperature field visualization
      spField = spVis->addScalarEquation<Equations::ScalarFieldVisualizer>();
      spField->setFields(true, false);
      spField->setIdentity(PhysicalNames::TEMPERATURE);

      // Add velocity fields visualization
      spField = spVis->addScalarEquation<Equations::ScalarFieldVisualizer>();
      spField->setFields(true, true);
      spField->setIdentity(PhysicalNames::VELOCITYX);
      spField = spVis->addScalarEquation<Equations::ScalarFieldVisualizer>();
      spField->setFields(true, true);
      spField->setIdentity(PhysicalNames::VELOCITYY);
      spField = spVis->addScalarEquation<Equations::ScalarFieldVisualizer>();
      spField->setFields(true, true);
      spField->setIdentity(PhysicalNames::VELOCITYZ);

      // Add output file
      IoVariable::SharedVisualizationFileWriter spOut(new IoVariable::VisualizationFileWriter(SchemeType::type()));
      spOut->expect(PhysicalNames::TEMPERATURE);
      spOut->expect(PhysicalNames::VELOCITY);
      spVis->addHdf5OutputFile(spOut);
   }

   void BoussinesqRB3DBoxVCModel::setVisualizationState(SharedVisualizationGenerator spVis)
   {
      // Create and add initial state file to IO
      IoVariable::SharedStateFileReader spIn(new IoVariable::StateFileReader("4Visu", SchemeType::type(), SchemeType::isRegular()));

      // Set expected fields
      spIn->expect(PhysicalNames::TEMPERATURE);
      spIn->expect(PhysicalNames::VELOCITY);

      // Set simulation state
      spVis->setInitialState(spIn);
   }

   void BoussinesqRB3DBoxVCModel::addAsciiOutputFiles(SharedSimulation spSim)
   {
      // Create maximal continuity writer
      IoVariable::SharedContinuityWriter spCont(new IoVariable::ContinuityWriter(SchemeType::type()));
      spCont->expect(PhysicalNames::VELOCITY);
      spSim->addAsciiOutputFile(spCont);

      // Create cubic Nusselt number writer
      IoVariable::SharedNusseltCubicWriter spNusselt(new IoVariable::NusseltCubicWriter(SchemeType::type()));
      spNusselt->expect(PhysicalNames::TEMPERATURE);
      spSim->addAsciiOutputFile(spNusselt);
   }

   void BoussinesqRB3DBoxVCModel::addHdf5OutputFiles(SharedSimulation spSim)
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

   void BoussinesqRB3DBoxVCModel::setInitialState(SharedSimulation spSim)
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
