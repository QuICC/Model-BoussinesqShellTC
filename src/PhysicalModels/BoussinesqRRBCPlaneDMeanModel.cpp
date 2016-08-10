/** 
 * @file BoussinesqRRBCPlaneDMeanModel.cpp
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
#include "PhysicalModels/BoussinesqRRBCPlaneDMeanModel.hpp"

// Project includes
//
#include "Enums/FieldIds.hpp"
#include "Equations/Box/Boussinesq/BoussinesqRRBCPlaneDMeanTransport.hpp"
#include "Equations/Box/Boussinesq/BoussinesqRRBCPlaneDMeanHeat.hpp"
#include "Equations/Box/Boussinesq/BoussinesqRRBCPlaneMomentum.hpp"
#include "IoVariable/StateFileReader.hpp"
#include "IoVariable/StateFileWriter.hpp"
#include "IoVariable/Cartesian1DScalarEnergyWriter.hpp"
#include "IoVariable/Cartesian1DTorPolEnergyWriter.hpp"
#include "IoVariable/Cartesian1DNusseltZWriter.hpp"
#include "IoVariable/VisualizationFileWriter.hpp"
#include "IoTools/IdToHuman.hpp"
#include "Generator/States/RandomScalarState.hpp"
#include "Generator/States/RandomVectorState.hpp"
#include "Generator/States/CartesianExactScalarState.hpp"
#include "Generator/States/CartesianExactVectorState.hpp"
#include "Generator/Visualizers/ScalarFieldVisualizer.hpp"
#include "Generator/Visualizers/VectorFieldVisualizer.hpp"
#include "PhysicalModels/PhysicalModelBase.hpp"

namespace GeoMHDiSCC {

   const std::string BoussinesqRRBCPlaneDMeanModel::PYMODULE = "boussinesq_rrbcplane_dmean";

   const std::string BoussinesqRRBCPlaneDMeanModel::PYCLASS = "BoussinesqRRBCPlane";

   const std::string BoussinesqRRBCPlaneDMeanModel::PYVISUCLASS = "BoussinesqRRBCPlaneVisu";

   void BoussinesqRRBCPlaneDMeanModel::addEquations(SharedSimulation spSim)
   {
      // Add transport equation
      spSim->addScalarEquation<Equations::BoussinesqRRBCPlaneDMeanTransport>();
      
      // Add Mean heat equation
      spSim->addScalarEquation<Equations::BoussinesqRRBCPlaneDMeanHeat>();
      
      // Add Navier-Stokes equation (X,Y,Z components)
      spSim->addVectorEquation<Equations::BoussinesqRRBCPlaneMomentum>();
   }

   void BoussinesqRRBCPlaneDMeanModel::addStates(SharedStateGenerator spGen)
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
         spVector->setStateType(Equations::CartesianExactStateIds::TORPOLTFF);

         // Add scalar exact initial state generator
         spScalar = spGen->addScalarEquation<Equations::CartesianExactScalarState>();
         spScalar->setIdentity(PhysicalNames::TEMPERATURE);
         spScalar->setStateType(Equations::CartesianExactStateIds::POLYSINSIN);
         spScalar->setModeOptions(1e0, 10.0, 1e0, 10.0, 1e0, 10.0);

      // Generate random spectrum
      } else
      {
         // Shared pointer to random initial state equation
         Equations::SharedRandomScalarState spScalar;
         Equations::SharedRandomVectorState spVector;

         // Add scalar random initial state generator 
         spVector = spGen->addVectorEquation<Equations::RandomVectorState>();
         spVector->setIdentity(PhysicalNames::VELOCITY);
         spVector->setSpectrum(FieldComponents::Spectral::TOR, -1e-6, 1e-6, 1e4, 1e4, 1e4);
         spVector->setSpectrum(FieldComponents::Spectral::POL, -1e-6, 1e-6, 1e4, 1e4, 1e4);

         // Add scalar random initial state generator
         spScalar = spGen->addScalarEquation<Equations::RandomScalarState>();
         spScalar->setIdentity(PhysicalNames::TEMPERATURE);
         spScalar->setSpectrum(-1e-6, 1e-6, 1e4, 1e4, 1e4, Equations::RandomScalarState::NOMEAN);

         // Add scalar random initial state generator
         spScalar = spGen->addScalarEquation<Equations::RandomScalarState>();
         spScalar->setIdentity(PhysicalNames::DZ_MEANTEMPERATURE);
         spScalar->setSpectrum(-1e-10, 1e-10, 1e4, 1e4, 1e4, Equations::RandomScalarState::ONLYMEAN);
      }

      // Add output file
      IoVariable::SharedStateFileWriter spOut(new IoVariable::StateFileWriter(SchemeType::type(), SchemeType::isRegular()));
      spOut->expect(PhysicalNames::VELOCITY);
      spOut->expect(PhysicalNames::TEMPERATURE);
      spOut->expect(PhysicalNames::DZ_MEANTEMPERATURE);
      spGen->addHdf5OutputFile(spOut);
   }

   void BoussinesqRRBCPlaneDMeanModel::addVisualizers(SharedVisualizationGenerator spVis)
   {
      // Shared pointer to basic field visualizer
      Equations::SharedScalarFieldVisualizer spScalar;
      Equations::SharedVectorFieldVisualizer spVector;

      // Add temperature field visualization
      spScalar = spVis->addScalarEquation<Equations::ScalarFieldVisualizer>();
      spScalar->setFields(true, false);
      spScalar->setIdentity(PhysicalNames::TEMPERATURE);

      // Add mean temperature gradient visualization
      spScalar = spVis->addScalarEquation<Equations::ScalarFieldVisualizer>();
      spScalar->setFields(true, false);
      spScalar->setIdentity(PhysicalNames::DZ_MEANTEMPERATURE);

      // Add velocity fields visualization
      spVector = spVis->addVectorEquation<Equations::VectorFieldVisualizer>();
      spVector->setFields(true, false, true);
      spVector->setIdentity(PhysicalNames::VELOCITY);

      // Add output file
      IoVariable::SharedVisualizationFileWriter spOut(new IoVariable::VisualizationFileWriter(SchemeType::type()));
      spOut->expect(PhysicalNames::TEMPERATURE);
      spOut->expect(PhysicalNames::DZ_MEANTEMPERATURE);
      spOut->expect(PhysicalNames::VELOCITY);
      spVis->addHdf5OutputFile(spOut);
   }

   void BoussinesqRRBCPlaneDMeanModel::setVisualizationState(SharedVisualizationGenerator spVis)
   {
      // Create and add initial state file to IO
      IoVariable::SharedStateFileReader spIn(new IoVariable::StateFileReader("4Visu", SchemeType::type(), SchemeType::isRegular()));

      // Set expected fields
      spIn->expect(PhysicalNames::TEMPERATURE);
      spIn->expect(PhysicalNames::DZ_MEANTEMPERATURE);
      spIn->expect(PhysicalNames::VELOCITY);

      // Set simulation state
      spVis->setInitialState(spIn);
   }

   void BoussinesqRRBCPlaneDMeanModel::addAsciiOutputFiles(SharedSimulation spSim)
   {
      // Create temperature energy writer
      IoVariable::SharedCartesian1DScalarEnergyWriter spTemp(new IoVariable::Cartesian1DScalarEnergyWriter("temperature", SchemeType::type()));
      spTemp->expect(PhysicalNames::TEMPERATURE);
      spSim->addAsciiOutputFile(spTemp);

      // Create temperature energy writer
      IoVariable::SharedCartesian1DScalarEnergyWriter spMeanTemp(new IoVariable::Cartesian1DScalarEnergyWriter("dz_meantemperature", SchemeType::type()));
      spMeanTemp->expect(PhysicalNames::DZ_MEANTEMPERATURE);
      spSim->addAsciiOutputFile(spMeanTemp);

      // Create kinetic energy writer
      IoVariable::SharedCartesian1DTorPolEnergyWriter spKinetic(new IoVariable::Cartesian1DTorPolEnergyWriter("kinetic", SchemeType::type()));
      spKinetic->expect(PhysicalNames::VELOCITY);
      spSim->addAsciiOutputFile(spKinetic);

      // Create nusselt number writer
      IoVariable::SharedCartesian1DNusseltZWriter spNusselt(new IoVariable::Cartesian1DNusseltZWriter(SchemeType::type()));
      spNusselt->expect(PhysicalNames::DZ_MEANTEMPERATURE);
      spSim->addAsciiOutputFile(spNusselt);
   }

   void BoussinesqRRBCPlaneDMeanModel::addHdf5OutputFiles(SharedSimulation spSim)
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

      spState->expect(PhysicalNames::DZ_MEANTEMPERATURE);

      spSim->addHdf5OutputFile(spState);
   }

   void BoussinesqRRBCPlaneDMeanModel::setInitialState(SharedSimulation spSim)
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
