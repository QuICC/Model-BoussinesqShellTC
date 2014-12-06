/** 
 * @file BoussinesqRRBAnnulusVCModel.cpp
 * @brief Source of the Boussinesq rotating Rayleigh-Benard annulus (velocity-continuity formulation) model
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
#include "PhysicalModels/BoussinesqRRBAnnulusVCModel.hpp"

// Project includes
//
#include "Enums/FieldIds.hpp"
#include "IoVariable/StateFileReader.hpp"
#include "IoVariable/StateFileWriter.hpp"
#include "IoVariable/ContinuityWriter.hpp"
#include "IoVariable/VisualizationFileWriter.hpp"
#include "IoTools/IdToHuman.hpp"
#include "Equations/Annulus/Boussinesq/BoussinesqRRBAnnulusVCTransport.hpp"
#include "Equations/Annulus/Boussinesq/BoussinesqRRBAnnulusVCVelocityX.hpp"
#include "Equations/Annulus/Boussinesq/BoussinesqRRBAnnulusVCVelocityY.hpp"
#include "Equations/Annulus/Boussinesq/BoussinesqRRBAnnulusVCVelocityZ.hpp"
#include "Equations/Annulus/Boussinesq/BoussinesqRRBAnnulusVCContinuity.hpp"
#include "Generator/States/RandomScalarState.hpp"
#include "Generator/States/AnnulusExactScalarState.hpp"
#include "Generator/Visualizers/ScalarFieldVisualizer.hpp"
#include "PhysicalModels/PhysicalModelBase.hpp"

namespace GeoMHDiSCC {

   const std::string BoussinesqRRBAnnulusVCModel::PYMODULE = "boussinesq_rrbannulus_vc";

   const std::string BoussinesqRRBAnnulusVCModel::PYCLASS = "BoussinesqRRBAnnulusVC";

   void BoussinesqRRBAnnulusVCModel::addEquations(SharedSimulation spSim)
   {
      // Add transport equation
      spSim->addScalarEquation<Equations::BoussinesqRRBAnnulusVCTransport>();
      
      // Add Navier-Stokes equation (X,Y,Z components)
      spSim->addScalarEquation<Equations::BoussinesqRRBAnnulusVCVelocityX>();
      spSim->addScalarEquation<Equations::BoussinesqRRBAnnulusVCVelocityY>();
      spSim->addScalarEquation<Equations::BoussinesqRRBAnnulusVCVelocityZ>();

      // Add continuity equation
      spSim->addScalarEquation<Equations::BoussinesqRRBAnnulusVCContinuity>();
   }

   void BoussinesqRRBAnnulusVCModel::addStates(SharedStateGenerator spGen)
   {
      // Generate "exact" solutions (trigonometric or monomial)
      if(false)
      {
         // Shared pointer to equation
         Equations::SharedAnnulusExactScalarState spExact;

         // Add scalar exact initial state generator
         spExact = spGen->addScalarEquation<Equations::AnnulusExactScalarState>();
         spExact->setIdentity(PhysicalNames::VELOCITYX);
         spExact->setStateType(Equations::AnnulusExactScalarState::POLYSINPOLY);
         spExact->setModeOptions(1e0, 0.0, 1e0, 1.0, 1e0, 1.0);

         // Add scalar exact initial state generator
         spExact = spGen->addScalarEquation<Equations::AnnulusExactScalarState>();
         spExact->setIdentity(PhysicalNames::VELOCITYY);
         spExact->setStateType(Equations::AnnulusExactScalarState::POLYCOSPOLY);
         spExact->setModeOptions(1e0, 1.0, 1e0, 1.0, 1e0, 1.0);

         // Add scalar exact initial state generator
         spExact = spGen->addScalarEquation<Equations::AnnulusExactScalarState>();
         spExact->setIdentity(PhysicalNames::VELOCITYZ);
         spExact->setStateType(Equations::AnnulusExactScalarState::POLYSINPOLY);
         spExact->setModeOptions(1e0, 1.0, 1e0, 3.0, 1e0, 0.0);

         // Add scalar exact initial state generator
         spExact = spGen->addScalarEquation<Equations::AnnulusExactScalarState>();
         spExact->setIdentity(PhysicalNames::TEMPERATURE);
         spExact->setStateType(Equations::AnnulusExactScalarState::POLYCOSPOLY);
         spExact->setModeOptions(1e0, 2.0, 1e0, 2.0, 1e0, 0.0);

      // Generate random spectrum
      } else
      {
         // Shared pointer to random initial state equation
         Equations::SharedRandomScalarState spRand;

         // Add scalar random initial state generator 
         spRand = spGen->addScalarEquation<Equations::RandomScalarState>();
         spRand->setIdentity(PhysicalNames::VELOCITYX);
         spRand->setSpectrum(-0.001, 0.001, 1e4, 1e4, 1e4);

         // Add scalar random initial state generator 
         spRand = spGen->addScalarEquation<Equations::RandomScalarState>();
         spRand->setIdentity(PhysicalNames::VELOCITYY);
         spRand->setSpectrum(-0.001, 0.001, 1e4, 1e4, 1e4);

         // Add scalar random initial state generator 
         spRand = spGen->addScalarEquation<Equations::RandomScalarState>();
         spRand->setIdentity(PhysicalNames::VELOCITYZ);
         spRand->setSpectrum(-0.001, 0.001, 1e4, 1e4, 1e4);

         // Add scalar random initial state generator
         spRand = spGen->addScalarEquation<Equations::RandomScalarState>();
         spRand->setIdentity(PhysicalNames::TEMPERATURE);
         spRand->setSpectrum(-0.001, 0.001, 1e4, 1e4, 1e4);
      }

      // Add output file
      IoVariable::SharedStateFileWriter spOut(new IoVariable::StateFileWriter(SchemeType::type(), SchemeType::isRegular()));
      spOut->expect(PhysicalNames::VELOCITYX);
      spOut->expect(PhysicalNames::VELOCITYY);
      spOut->expect(PhysicalNames::VELOCITYZ);
      spOut->expect(PhysicalNames::TEMPERATURE);
      spGen->addHdf5OutputFile(spOut);
   }

   void BoussinesqRRBAnnulusVCModel::addVisualizers(SharedVisualizationGenerator spVis)
   {
      // Shared pointer to basic field visualizer
      Equations::SharedScalarFieldVisualizer spField;

      // Add temperature field visualization
      spField = spVis->addScalarEquation<Equations::ScalarFieldVisualizer>();
      spField->setFields(true, true);
      spField->setIdentity(PhysicalNames::TEMPERATURE);

      // Add velocity field visualization
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
      spOut->expect(PhysicalNames::VELOCITYX);
      spOut->expect(PhysicalNames::VELOCITYY);
      spOut->expect(PhysicalNames::VELOCITYZ);
      spVis->addHdf5OutputFile(spOut);
   }

   void BoussinesqRRBAnnulusVCModel::setVisualizationState(SharedVisualizationGenerator spVis)
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

   void BoussinesqRRBAnnulusVCModel::addAsciiOutputFiles(SharedSimulation spSim)
   {
      // Create maximal continuity writer
      IoVariable::SharedContinuityWriter spState(new IoVariable::ContinuityWriter(SchemeType::type()));
      spState->expect(PhysicalNames::VELOCITYX);
      spState->expect(PhysicalNames::VELOCITYY);
      spState->expect(PhysicalNames::VELOCITYZ);
      spSim->addAsciiOutputFile(spState);
   }

   void BoussinesqRRBAnnulusVCModel::addHdf5OutputFiles(SharedSimulation spSim)
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

   void BoussinesqRRBAnnulusVCModel::setInitialState(SharedSimulation spSim)
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
