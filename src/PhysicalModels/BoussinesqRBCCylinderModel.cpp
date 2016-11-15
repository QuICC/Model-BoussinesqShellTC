/** 
 * @file BoussinesqRBCCylinderModel.cpp
 * @brief Source of Boussinesq Rayleigh-Benard convection in a cylinder (toroidal-poloidal formulation)
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
#include "PhysicalModels/BoussinesqRBCCylinderModel.hpp"

// Project includes
//
#include "Enums/FieldIds.hpp"
#include "Equations/Cylinder/Boussinesq/BoussinesqRBCCylinderTransport.hpp"
#include "Equations/Cylinder/Boussinesq/BoussinesqRBCCylinderMomentum.hpp"
#include "IoVariable/StateFileReader.hpp"
#include "IoVariable/StateFileWriter.hpp"
#include "IoVariable/VisualizationFileWriter.hpp"
#include "IoTools/IdToHuman.hpp"
#include "IoVariable/CylinderScalarEnergyWriter.hpp"
#include "IoVariable/CylinderTorPolEnergyWriter.hpp"
#include "Generator/States/RandomScalarState.hpp"
#include "Generator/States/RandomVectorState.hpp"
#include "Generator/States/CylinderExactStateIds.hpp"
#include "Generator/States/CylinderExactScalarState.hpp"
#include "Generator/States/CylinderExactVectorState.hpp"
#include "Generator/Visualizers/ScalarFieldVisualizer.hpp"
#include "Generator/Visualizers/ScalarFieldTrivialVisualizer.hpp"
#include "Generator/Visualizers/VectorFieldVisualizer.hpp"
#include "Generator/Visualizers/VectorFieldTrivialVisualizer.hpp"
#include "PhysicalModels/PhysicalModelBase.hpp"

namespace GeoMHDiSCC {

   const std::string BoussinesqRBCCylinderModel::PYMODULE = "boussinesq_rbccylinder";

   const std::string BoussinesqRBCCylinderModel::PYCLASS = "BoussinesqRBCCylinder";

   void BoussinesqRBCCylinderModel::addEquations(SharedSimulation spSim)
   {
      // Add transport equation
      spSim->addScalarEquation<Equations::BoussinesqRBCCylinderTransport>();
      
      // Add Navier-Stokes equation (X,Y,Z components)
      spSim->addVectorEquation<Equations::BoussinesqRBCCylinderMomentum>();
   }

   void BoussinesqRBCCylinderModel::addStates(SharedStateGenerator spGen)
   {
      // Generate "exact" solutions (trigonometric or monomial)
      if(false)
      {
         // Shared pointer to equation
         Equations::SharedCylinderExactScalarState spScalar;
         Equations::SharedCylinderExactVectorState spVector;

         // Add scalar exact initial state generator
         spVector = spGen->addVectorEquation<Equations::CylinderExactVectorState>();
         spVector->setIdentity(PhysicalNames::VELOCITY);
         spVector->setStateType(FieldComponents::Physical::R, Equations::CylinderExactStateIds::POLYCOSPOLY);
         spVector->setModeOptions(FieldComponents::Physical::R, 1.0e0, 1.0, 1.0e0, 0.0, 1.0e0, 0.0);
         spVector->setStateType(FieldComponents::Physical::THETA, Equations::CylinderExactStateIds::POLYCOSPOLY);
         spVector->setModeOptions(FieldComponents::Physical::THETA, 1.0e0, 0.0, 1.0e0, 1.0, 1.0e0, 0.0);
         spVector->setStateType(FieldComponents::Physical::Z, Equations::CylinderExactStateIds::POLYCOSPOLY);
         spVector->setModeOptions(FieldComponents::Physical::Z, 1.0e0, 0.0, 1.0e0, 0.0, 1.0e0, 1.0);

         // Add scalar exact initial state generator
         spScalar = spGen->addScalarEquation<Equations::CylinderExactScalarState>();
         spScalar->setIdentity(PhysicalNames::TEMPERATURE);
         spScalar->setStateType(Equations::CylinderExactStateIds::POLYCOSPOLY);
         spScalar->setModeOptions(1e0, 2.0, 1e0, 2.0, 1e0, 1.0);

      // Generate random spectrum
      } else
      {
         // Shared pointer to random initial state equation
         Equations::SharedRandomScalarState spScalar;
         Equations::SharedRandomVectorState spVector;

         // Add scalar random initial state generator 
         spVector = spGen->addVectorEquation<Equations::RandomVectorState>();
         spVector->setIdentity(PhysicalNames::VELOCITY);
         spVector->setSpectrum(FieldComponents::Spectral::TOR, -1e-10, 1e-10, 1e4, 1e4, 1e4);
         spVector->setSpectrum(FieldComponents::Spectral::POL, -1e-10, 1e-10, 1e4, 1e4, 1e4);

         // Add scalar random initial state generator
         spScalar = spGen->addScalarEquation<Equations::RandomScalarState>();
         spScalar->setIdentity(PhysicalNames::TEMPERATURE);
         spScalar->setSpectrum(-1e-10, 1e-10, 1e4, 1e4, 1e4);
      }

      // Add output file
      IoVariable::SharedStateFileWriter spOut(new IoVariable::StateFileWriter(SchemeType::type(), SchemeType::isRegular()));
      spOut->expect(PhysicalNames::VELOCITY);
      spOut->expect(PhysicalNames::TEMPERATURE);
      spGen->addHdf5OutputFile(spOut);
   }

   void BoussinesqRBCCylinderModel::addVisualizers(SharedVisualizationGenerator spVis)
   {
      // Shared pointer to basic field visualizer
      Equations::SharedScalarFieldVisualizer spScalar;
      Equations::SharedScalarFieldTrivialVisualizer spSTrivial;
      Equations::SharedVectorFieldVisualizer spVector;
      Equations::SharedVectorFieldTrivialVisualizer spVTrivial;

      // Add temperature field visualization
      spScalar = spVis->addScalarEquation<Equations::ScalarFieldVisualizer>();
      spScalar->setFields(true, false);
      spScalar->setIdentity(PhysicalNames::TEMPERATURE);

      // Add mean temperature field visualization
      spSTrivial = spVis->addScalarEquation<Equations::ScalarFieldTrivialVisualizer>();
      spSTrivial->setFields(true, false);
      spSTrivial->setIdentity(PhysicalNames::MEAN_TEMPERATURE);

      // Add fluctuating temperature field visualization
      spSTrivial = spVis->addScalarEquation<Equations::ScalarFieldTrivialVisualizer>();
      spSTrivial->setFields(true, false);
      spSTrivial->setIdentity(PhysicalNames::FLUCT_TEMPERATURE);

      // Add velocity fields visualization
      spVector = spVis->addVectorEquation<Equations::VectorFieldVisualizer>();
      spVector->setFields(true, false, true);
      spVector->setIdentity(PhysicalNames::VELOCITY);

      // Add velocity fields visualization
      spVTrivial = spVis->addVectorEquation<Equations::VectorFieldTrivialVisualizer>();
      spVTrivial->setFields(true, false, true);
      spVTrivial->setIdentity(PhysicalNames::MAGNETIC);

      // Add output file
      IoVariable::SharedVisualizationFileWriter spOut(new IoVariable::VisualizationFileWriter(SchemeType::type()));
      spOut->expect(PhysicalNames::TEMPERATURE);
      spOut->expect(PhysicalNames::MEAN_TEMPERATURE);
      spOut->expect(PhysicalNames::FLUCT_TEMPERATURE);
      spOut->expect(PhysicalNames::VELOCITY);
      spOut->expect(PhysicalNames::MAGNETIC);
      spVis->addHdf5OutputFile(spOut);
   }

   void BoussinesqRBCCylinderModel::setVisualizationState(SharedVisualizationGenerator spVis)
   {
      // Create and add initial state file to IO
      IoVariable::SharedStateFileReader spIn(new IoVariable::StateFileReader("4Visu", SchemeType::type(), SchemeType::isRegular()));

      // Set expected fields
      spIn->expect(PhysicalNames::TEMPERATURE);
      spIn->expect(PhysicalNames::VELOCITY);

      // Set simulation state
      spVis->setInitialState(spIn);
   }

   void BoussinesqRBCCylinderModel::addAsciiOutputFiles(SharedSimulation spSim)
   {
      // Create temperature energy writer
      IoVariable::SharedCylinderScalarEnergyWriter spScalar(new IoVariable::CylinderScalarEnergyWriter("temperature", SchemeType::type()));
      spScalar->expect(PhysicalNames::TEMPERATURE);
      spSim->addAsciiOutputFile(spScalar);

      // Create kinetic energy writer
      IoVariable::SharedCylinderTorPolEnergyWriter spVector(new IoVariable::CylinderTorPolEnergyWriter("kinetic", SchemeType::type()));
      spVector->expect(PhysicalNames::VELOCITY);
      spSim->addAsciiOutputFile(spVector);
   }

   void BoussinesqRBCCylinderModel::addHdf5OutputFiles(SharedSimulation spSim)
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

   void BoussinesqRBCCylinderModel::setInitialState(SharedSimulation spSim)
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
