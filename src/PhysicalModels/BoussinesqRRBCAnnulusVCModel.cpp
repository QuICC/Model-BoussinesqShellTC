/** 
 * @file BoussinesqRRBCAnnulusVCModel.cpp
 * @brief Source of the Boussinesq rotating Rayleigh-Benard convection annulus (velocity-continuity formulation) model
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
#include "PhysicalModels/BoussinesqRRBCAnnulusVCModel.hpp"

// Project includes
//
#include "Enums/FieldIds.hpp"
#include "IoVariable/StateFileReader.hpp"
#include "IoVariable/StateFileWriter.hpp"
#include "IoVariable/ContinuityWriter.hpp"
#include "IoVariable/VisualizationFileWriter.hpp"
#include "IoTools/IdToHuman.hpp"
#include "Equations/Annulus/Boussinesq/BoussinesqRRBCAnnulusVCTransport.hpp"
#include "Equations/Annulus/Boussinesq/BoussinesqRRBCAnnulusVCMomentum.hpp"
#include "Equations/Annulus/Boussinesq/BoussinesqRRBCAnnulusVCContinuity.hpp"
#include "Generator/States/RandomScalarState.hpp"
#include "Generator/States/RandomVectorState.hpp"
#include "Generator/States/AnnulusExactScalarState.hpp"
#include "Generator/States/AnnulusExactVectorState.hpp"
#include "Generator/Visualizers/ScalarFieldVisualizer.hpp"
#include "Generator/Visualizers/VectorFieldVisualizer.hpp"
#include "PhysicalModels/PhysicalModelBase.hpp"

namespace QuICC {

   const std::string BoussinesqRRBCAnnulusVCModel::PYMODULE = "boussinesq_rrbcannulus_vc";

   const std::string BoussinesqRRBCAnnulusVCModel::PYCLASS = "BoussinesqRRBCAnnulusVC";

   void BoussinesqRRBCAnnulusVCModel::addEquations(SharedSimulation spSim)
   {
      // Add transport equation
      spSim->addScalarEquation<Equations::BoussinesqRRBCAnnulusVCTransport>();
      
      // Add Navier-Stokes equation (X,Y,Z components)
      spSim->addVectorEquation<Equations::BoussinesqRRBCAnnulusVCMomentum>();

      // Add continuity equation
      spSim->addScalarEquation<Equations::BoussinesqRRBCAnnulusVCContinuity>();
   }

   void BoussinesqRRBCAnnulusVCModel::addStates(SharedStateGenerator spGen)
   {
      // Generate "exact" solutions (trigonometric or monomial)
      if(true)
      {
         // Shared pointer to equation
         Equations::SharedAnnulusExactScalarState spScalar;
         Equations::SharedAnnulusExactVectorState spVector;

         // Add scalar exact initial state generator
         spVector = spGen->addVectorEquation<Equations::AnnulusExactVectorState>();
         spVector->setIdentity(PhysicalNames::VELOCITY);
         spVector->setStateType(FieldComponents::Physical::R, Equations::AnnulusExactStateIds::POLYCOSPOLY);
         spVector->setModeOptions(FieldComponents::Physical::R, 1.0e0, 1.0, 1.0e0, 0.0, 1.0e0, 0.0);
         spVector->setStateType(FieldComponents::Physical::THETA, Equations::AnnulusExactStateIds::POLYCOSPOLY);
         spVector->setModeOptions(FieldComponents::Physical::THETA, 1.0e0, 0.0, 1.0e0, 1.0, 1.0e0, 0.0);
         spVector->setStateType(FieldComponents::Physical::Z, Equations::AnnulusExactStateIds::POLYCOSPOLY);
         spVector->setModeOptions(FieldComponents::Physical::Z, 1.0e0, 0.0, 1.0e0, 0.0, 1.0e0, 1.0);

         // Add scalar exact initial state generator
         spScalar = spGen->addScalarEquation<Equations::AnnulusExactScalarState>();
         spScalar->setIdentity(PhysicalNames::TEMPERATURE);
         spScalar->setStateType(Equations::AnnulusExactStateIds::POLYCOSPOLY);
         spScalar->setModeOptions(1e0, 2.0, 1e0, 2.0, 1e0, 0.0);

      // Generate random spectrum
      } else
      {
         // Shared pointer to random initial state equation
         Equations::SharedRandomScalarState spScalar;
         Equations::SharedRandomVectorState spVector;

         // Add scalar random initial state generator 
         spVector = spGen->addVectorEquation<Equations::RandomVectorState>();
         spVector->setIdentity(PhysicalNames::VELOCITY);
         spVector->setSpectrum(FieldComponents::Spectral::R, -1e-4, 1e-4, 1e4, 1e4, 1e4);
         spVector->setSpectrum(FieldComponents::Spectral::THETA, -1e-4, 1e-4, 1e4, 1e4, 1e4);
         spVector->setSpectrum(FieldComponents::Spectral::Z, -1e-4, 1e-4, 1e4, 1e4, 1e4);

         // Add scalar random initial state generator
         spScalar = spGen->addScalarEquation<Equations::RandomScalarState>();
         spScalar->setIdentity(PhysicalNames::TEMPERATURE);
         spScalar->setSpectrum(-1e-3, 1e3, 1e4, 1e4, 1e4);
      }

      // Add output file
      IoVariable::SharedStateFileWriter spOut(new IoVariable::StateFileWriter(SchemeType::type(), SchemeType::isRegular()));
      spOut->expect(PhysicalNames::VELOCITY);
      spOut->expect(PhysicalNames::TEMPERATURE);
      spGen->addHdf5OutputFile(spOut);
   }

   void BoussinesqRRBCAnnulusVCModel::addVisualizers(SharedVisualizationGenerator spVis)
   {
      // Shared pointer to basic field visualizer
      Equations::SharedScalarFieldVisualizer spScalar;
      Equations::SharedVectorFieldVisualizer spVector;

      // Add temperature field visualization
      spScalar = spVis->addScalarEquation<Equations::ScalarFieldVisualizer>();
      spScalar->setFields(true, true);
      spScalar->setIdentity(PhysicalNames::TEMPERATURE);

      // Add velocity fields visualization
      spVector = spVis->addVectorEquation<Equations::VectorFieldVisualizer>();
      spVector->setFields(true, true, false);
      spVector->setIdentity(PhysicalNames::VELOCITY);

      // Add output file
      IoVariable::SharedVisualizationFileWriter spOut(new IoVariable::VisualizationFileWriter(SchemeType::type()));
      spOut->expect(PhysicalNames::TEMPERATURE);
      spOut->expect(PhysicalNames::VELOCITY);
      spVis->addHdf5OutputFile(spOut);
   }

   void BoussinesqRRBCAnnulusVCModel::setVisualizationState(SharedVisualizationGenerator spVis)
   {
      // Create and add initial state file to IO
      IoVariable::SharedStateFileReader spIn(new IoVariable::StateFileReader("4Visu", SchemeType::type(), SchemeType::isRegular()));

      // Set expected fields
      spIn->expect(PhysicalNames::TEMPERATURE);
      spIn->expect(PhysicalNames::VELOCITY);

      // Set simulation state
      spVis->setInitialState(spIn);
   }

   void BoussinesqRRBCAnnulusVCModel::addAsciiOutputFiles(SharedSimulation spSim)
   {
      // Create maximal continuity writer
//      IoVariable::SharedContinuityWriter spState(new IoVariable::ContinuityWriter(SchemeType::type()));
//      spState->expect(PhysicalNames::VELOCITY);
//      spSim->addAsciiOutputFile(spState);
   }

   void BoussinesqRRBCAnnulusVCModel::addHdf5OutputFiles(SharedSimulation spSim)
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
      spSim->addHdf5OutputFile(spState);
   }

   void BoussinesqRRBCAnnulusVCModel::addStatsOutputFiles(SharedSimulation spSim)
   {
   }

   void BoussinesqRRBCAnnulusVCModel::setInitialState(SharedSimulation spSim)
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
