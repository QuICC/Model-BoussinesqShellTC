/** 
 * @file BoussinesqRRBCPlaneVCModel.cpp
 * @brief Source of the rotating Boussinesq Rayleigh-Benard convection in a plane layer (velocity-continuity formulation) model
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
#include "PhysicalModels/BoussinesqRRBCPlaneVCModel.hpp"

// Project includes
//
#include "Enums/FieldIds.hpp"
#include "Equations/Box/Boussinesq/BoussinesqRRBCPlaneVCTransport.hpp"
#include "Equations/Box/Boussinesq/BoussinesqRRBCPlaneVCMomentum.hpp"
#include "Equations/Box/Boussinesq/BoussinesqRRBCPlaneVCContinuity.hpp"
#include "IoVariable/StateFileReader.hpp"
#include "IoVariable/StateFileWriter.hpp"
#include "IoVariable/ContinuityWriter.hpp"
#include "IoVariable/Cartesian1DScalarEnergyWriter.hpp"
#include "IoVariable/Cartesian1DPrimitiveEnergyWriter.hpp"
#include "IoVariable/VisualizationFileWriter.hpp"
#include "IoTools/IdToHuman.hpp"
#include "Generator/States/RandomScalarState.hpp"
#include "Generator/States/RandomVectorState.hpp"
#include "Generator/States/CartesianExactScalarState.hpp"
#include "Generator/States/CartesianExactVectorState.hpp"
#include "Generator/Visualizers/ScalarFieldVisualizer.hpp"
#include "Generator/Visualizers/VectorFieldVisualizer.hpp"
#include "PhysicalModels/PhysicalModelBase.hpp"

namespace QuICC {

   const std::string BoussinesqRRBCPlaneVCModel::PYMODULE = "boussinesq_rrbcplane_vc";

   const std::string BoussinesqRRBCPlaneVCModel::PYCLASS = "BoussinesqRRBCPlaneVC";

   void BoussinesqRRBCPlaneVCModel::addEquations(SharedSimulation spSim)
   {
      // Add transport equation
      spSim->addScalarEquation<Equations::BoussinesqRRBCPlaneVCTransport>();
      
      // Add Navier-Stokes equation (X,Y,Z components)
      spSim->addVectorEquation<Equations::BoussinesqRRBCPlaneVCMomentum>();

      // Add continuity equation
      spSim->addScalarEquation<Equations::BoussinesqRRBCPlaneVCContinuity>();
   }

   void BoussinesqRRBCPlaneVCModel::addStates(SharedStateGenerator spGen)
   {
      // Generate "exact" solutions (trigonometric or monomial)
      if(false)
      {
         // Shared pointer to equation
         Equations::SharedCartesianExactScalarState spScalar;
         Equations::SharedCartesianExactVectorState spVector;

         // Add vector exact initial state generator
         spVector = spGen->addVectorEquation<Equations::CartesianExactVectorState>();
         spVector->setIdentity(PhysicalNames::VELOCITY);
         spVector->setStateType(FieldComponents::Physical::X, Equations::CartesianExactStateIds::POLYCOSCOS);
         spVector->setModeOptions(FieldComponents::Physical::X, 1.0e0, 1.0, 1.0e0, 1.0, 1.0e0, 1.0);
         spVector->setStateType(FieldComponents::Physical::Y, Equations::CartesianExactStateIds::POLYCOSCOS);
         spVector->setModeOptions(FieldComponents::Physical::Y, 1.0e0, 1.0, 1.0e0, 1.0, 1.0e0, 1.0);
         spVector->setStateType(FieldComponents::Physical::Z, Equations::CartesianExactStateIds::POLYCOSCOS);
         spVector->setModeOptions(FieldComponents::Physical::Z, 1.0e0, 1.0, 1.0e0, 1.0, 1.0e0, 1.0);

         // Add scalar exact initial state generator
         spScalar = spGen->addScalarEquation<Equations::CartesianExactScalarState>();
         spScalar->setIdentity(PhysicalNames::TEMPERATURE);
         spScalar->setStateType(Equations::CartesianExactStateIds::POLYSINSIN);
         spScalar->setModeOptions(-1e2, 10.0, 3e0, 10.0, -3e1, 10.0);

      // Generate random spectrum
      } else
      {
         // Shared pointer to random initial state equation
         Equations::SharedRandomScalarState spScalar;
         Equations::SharedRandomVectorState spVector;

         // Add scalar random initial state generator 
         spVector = spGen->addVectorEquation<Equations::RandomVectorState>();
         spVector->setIdentity(PhysicalNames::VELOCITY);
         spVector->setSpectrum(FieldComponents::Spectral::X, -1e-4, 1e-4, 1e4, 1e4, 1e4);
         spVector->setSpectrum(FieldComponents::Spectral::Y, -1e-4, 1e-4, 1e4, 1e4, 1e4);
         spVector->setSpectrum(FieldComponents::Spectral::Z, -1e-4, 1e-4, 1e4, 1e4, 1e4);

         // Add scalar random initial state generator
         spScalar = spGen->addScalarEquation<Equations::RandomScalarState>();
         spScalar->setIdentity(PhysicalNames::TEMPERATURE);
         spScalar->setSpectrum(-1e-4, 1e-4, 1e4, 1e4, 1e4);
      }

      // Add output file
      IoVariable::SharedStateFileWriter spOut(new IoVariable::StateFileWriter(SchemeType::type(), SchemeType::isRegular()));
      spOut->expect(PhysicalNames::VELOCITY);
      spOut->expect(PhysicalNames::TEMPERATURE);
      spGen->addHdf5OutputFile(spOut);
   }

   void BoussinesqRRBCPlaneVCModel::addVisualizers(SharedVisualizationGenerator spVis)
   {
      // Shared pointer to basic field visualizer
      Equations::SharedScalarFieldVisualizer spScalar;
      Equations::SharedVectorFieldVisualizer spVector;

      // Add temperature field visualization
      spScalar = spVis->addScalarEquation<Equations::ScalarFieldVisualizer>();
      spScalar->setFields(true, false);
      spScalar->setIdentity(PhysicalNames::TEMPERATURE);

      // Add velocity fields visualization
      spVector = spVis->addVectorEquation<Equations::VectorFieldVisualizer>();
      spVector->setFields(true, false, false);
      spVector->setIdentity(PhysicalNames::VELOCITY);

      // Add pressure field visualization
      spScalar = spVis->addScalarEquation<Equations::ScalarFieldVisualizer>();
      spScalar->setFields(true, false);
      spScalar->setIdentity(PhysicalNames::PRESSURE);

      // Add output file
      IoVariable::SharedVisualizationFileWriter spOut(new IoVariable::VisualizationFileWriter(SchemeType::type()));
      spOut->expect(PhysicalNames::TEMPERATURE);
      spOut->expect(PhysicalNames::VELOCITY);
      spOut->expect(PhysicalNames::PRESSURE);
      spVis->addHdf5OutputFile(spOut);
   }

   void BoussinesqRRBCPlaneVCModel::setVisualizationState(SharedVisualizationGenerator spVis)
   {
      // Create and add initial state file to IO
      IoVariable::SharedStateFileReader spIn(new IoVariable::StateFileReader("4Visu", SchemeType::type(), SchemeType::isRegular()));

      // Set expected fields
      spIn->expect(PhysicalNames::TEMPERATURE);
      spIn->expect(PhysicalNames::VELOCITY);
      spIn->expect(PhysicalNames::PRESSURE, false);

      // Set simulation state
      spVis->setInitialState(spIn);
   }

   void BoussinesqRRBCPlaneVCModel::addAsciiOutputFiles(SharedSimulation spSim)
   {
      // Create maximal continuity writer
      IoVariable::SharedContinuityWriter spCont(new IoVariable::ContinuityWriter(SchemeType::type()));
      spCont->expect(PhysicalNames::VELOCITY);
      spSim->addAsciiOutputFile(spCont);

      // Create temperature energy writer
      IoVariable::SharedCartesian1DScalarEnergyWriter spTemp(new IoVariable::Cartesian1DScalarEnergyWriter("temperature", SchemeType::type()));
      spTemp->expect(PhysicalNames::TEMPERATURE);
      spSim->addAsciiOutputFile(spTemp);

      // Create kinetic energy writer
      IoVariable::SharedCartesian1DPrimitiveEnergyWriter spKinetic(new IoVariable::Cartesian1DPrimitiveEnergyWriter("kinetic", SchemeType::type()));
      spKinetic->expect(PhysicalNames::VELOCITY);
      spSim->addAsciiOutputFile(spKinetic);
   }

   void BoussinesqRRBCPlaneVCModel::addHdf5OutputFiles(SharedSimulation spSim)
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
      spState->expect(PhysicalNames::PRESSURE);
      spSim->addHdf5OutputFile(spState);
   }

   void BoussinesqRRBCPlaneVCModel::addStatsOutputFiles(SharedSimulation spSim)
   {
   }

   void BoussinesqRRBCPlaneVCModel::setInitialState(SharedSimulation spSim)
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