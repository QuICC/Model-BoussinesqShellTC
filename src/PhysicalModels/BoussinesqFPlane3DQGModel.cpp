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
#include "Equations/Asymptotics/FPlane3DQG/Boussinesq/BoussinesqFPlane3DQGStreamfunction.hpp"
#include "Equations/Asymptotics/FPlane3DQG/Boussinesq/BoussinesqFPlane3DQGVelocityZ.hpp"
#include "Equations/Asymptotics/FPlane3DQG/Boussinesq/BoussinesqFPlane3DQGTransport.hpp"
#include "Equations/Asymptotics/FPlane3DQG/Boussinesq/BoussinesqFPlane3DQGVorticityZ.hpp"
#include "Equations/Asymptotics/FPlane3DQG/Boussinesq/BoussinesqFPlane3DQGMeanHeat.hpp"
#include "IoVariable/Cartesian1DNusseltZWriter.hpp"
#include "IoVariable/Cartesian1DScalarEnergyWriter.hpp"
#include "IoVariable/Cartesian1DStreamEnergyWriter.hpp"
#include "Generator/States/RandomScalarState.hpp"
#include "Generator/States/CartesianExactScalarState.hpp"
#include "Generator/Visualizers/ScalarFieldVisualizer.hpp"
#include "PhysicalModels/PhysicalModelBase.hpp"

namespace GeoMHDiSCC {

   const std::string BoussinesqFPlane3DQGModel::PYMODULE = "boussinesq_fplane3dqg";

   const std::string BoussinesqFPlane3DQGModel::PYCLASS = "BoussinesqFPlane3DQG";

   void BoussinesqFPlane3DQGModel::addEquations(SharedSimulation spSim)
   {
      // Add upright streamfunction equation
      spSim->addScalarEquation<Equations::BoussinesqFPlane3DQGStreamfunction>();
      
      // Add upright vertical velocity equation
      spSim->addScalarEquation<Equations::BoussinesqFPlane3DQGVelocityZ>();
      
      // Add upright transport equation
      spSim->addScalarEquation<Equations::BoussinesqFPlane3DQGTransport>();


      // Add vertical vorticity equation
      spSim->addScalarEquation<Equations::BoussinesqFPlane3DQGVorticityZ>();

      
      // Add mean heat computation
      spSim->addScalarEquation<Equations::BoussinesqFPlane3DQGMeanHeat>();
   }

   void BoussinesqFPlane3DQGModel::addStates(SharedStateGenerator spGen)
   {
      // Generate "exact" solutions (trigonometric or monomial)
      if(false)
      {
         // Shared pointer to equation
         Equations::SharedCartesianExactScalarState spExact;

         // Add transport initial state generation equation
         spExact = spGen->addScalarEquation<Equations::CartesianExactScalarState>();
         spExact->setIdentity(PhysicalNames::TEMPERATURE);
         spExact->setStateType(Equations::CartesianExactStateIds::POLYCOSCOS);
         spExact->setModeOptions(1e0, 0.0, 1e0, 0.0, 1e0, 0.0);

         // Add streamfunction initial state generation equation
         spExact = spGen->addScalarEquation<Equations::CartesianExactScalarState>();
         spExact->setIdentity(PhysicalNames::STREAMFUNCTION);
         spExact->setStateType(Equations::CartesianExactStateIds::POLYCOSCOS);
         spExact->setModeOptions(2e0, 0.0, 1e0, 0.0, 1e0, 0.0);

         // Add vertical velocity initial state generation equation
         spExact = spGen->addScalarEquation<Equations::CartesianExactScalarState>();
         spExact->setIdentity(PhysicalNames::VELOCITYZ);
         spExact->setStateType(Equations::CartesianExactStateIds::POLYCOSCOS);
         spExact->setModeOptions(3e0, 0.0, 1e0, 0.0, 1e0, 0.0);

         // Add vertical velocity initial state generation equation
         spExact = spGen->addScalarEquation<Equations::CartesianExactScalarState>();
         spExact->setIdentity(PhysicalNames::DZ_MEANTEMPERATURE);
         spExact->setStateType(Equations::CartesianExactStateIds::POLYCOSCOS);
         spExact->setModeOptions(-1e0, 0.0, 1e0, 0.0, 1e0, 0.0);

      // Generate random spectrum
      } else
      {
         // Shared pointer to random initial state equation
         Equations::SharedRandomScalarState spRand;
         Equations::SharedCartesianExactScalarState spExact;

         // Add transport initial state generation equation
         spRand = spGen->addScalarEquation<Equations::RandomScalarState>();
         spRand->setIdentity(PhysicalNames::TEMPERATURE);
         spRand->setSpectrum(-1e-2, 1e-2, 1e4, 1e4, 1e4);

         // Add streamfunction initial state generation equation
         spRand = spGen->addScalarEquation<Equations::RandomScalarState>();
         spRand->setIdentity(PhysicalNames::STREAMFUNCTION);
         spRand->setSpectrum(-1e-2, 1e-2, 1e4, 1e4, 1e4);

         // Add vertical velocity initial state generation equation
         spRand = spGen->addScalarEquation<Equations::RandomScalarState>();
         spRand->setIdentity(PhysicalNames::VELOCITYZ);
         spRand->setSpectrum(-1e-2, 1e-2, 1e4, 1e4, 1e4);

         // Add vertical velocity initial state generation equation
         spExact = spGen->addScalarEquation<Equations::CartesianExactScalarState>();
         spExact->setIdentity(PhysicalNames::DZ_MEANTEMPERATURE);
         spExact->setStateType(Equations::CartesianExactStateIds::POLYCOSCOS);
         spExact->setModeOptions(0e0, 0.0, 1e0, 0.0, 1e0, 0.0);
      }

      // Add output file
      IoVariable::SharedStateFileWriter spOut(new IoVariable::StateFileWriter(SchemeType::type(), SchemeType::isRegular()));
      spOut->expect(PhysicalNames::TEMPERATURE);
      spOut->expect(PhysicalNames::STREAMFUNCTION);
      spOut->expect(PhysicalNames::VELOCITYZ);
      spOut->expect(PhysicalNames::DZ_MEANTEMPERATURE);
      spGen->addHdf5OutputFile(spOut);
   }

   void BoussinesqFPlane3DQGModel::addVisualizers(SharedVisualizationGenerator spVis)
   {
      // Shared pointer to basic field visualizer
      Equations::SharedScalarFieldVisualizer spField;

      // Add transport field visualization
      spField = spVis->addScalarEquation<Equations::ScalarFieldVisualizer>();
      spField->setFields(true, false);
      spField->setIdentity(PhysicalNames::TEMPERATURE);
      
      // Add streamfunction field visualization
      spField = spVis->addScalarEquation<Equations::ScalarFieldVisualizer>();
      spField->setFields(true, false);
      spField->setIdentity(PhysicalNames::STREAMFUNCTION);
      
      // Add vertical velocity field visualization
      spField = spVis->addScalarEquation<Equations::ScalarFieldVisualizer>();
      spField->setFields(true, false);
      spField->setIdentity(PhysicalNames::VELOCITYZ);
      
      // Add background temperature profile visualization
      spField = spVis->addScalarEquation<Equations::ScalarFieldVisualizer>();
      spField->setFields(true, false);
      spField->setIdentity(PhysicalNames::DZ_MEANTEMPERATURE);

      // Add output file
      IoVariable::SharedVisualizationFileWriter spOut(new IoVariable::VisualizationFileWriter(SchemeType::type()));
      spOut->expect(PhysicalNames::TEMPERATURE);
      spOut->expect(PhysicalNames::STREAMFUNCTION);
      spOut->expect(PhysicalNames::VELOCITYZ);
      spOut->expect(PhysicalNames::DZ_MEANTEMPERATURE);
      spVis->addHdf5OutputFile(spOut);
   }

   void BoussinesqFPlane3DQGModel::setVisualizationState(SharedVisualizationGenerator spVis)
   {
      // Create and add initial state file to IO
      IoVariable::SharedStateFileReader spIn(new IoVariable::StateFileReader("4Visu", SchemeType::type(), SchemeType::isRegular()));

      // Set expected fields
      spIn->expect(PhysicalNames::TEMPERATURE);
      spIn->expect(PhysicalNames::STREAMFUNCTION);
      spIn->expect(PhysicalNames::VELOCITYZ);
      spIn->expect(PhysicalNames::DZ_MEANTEMPERATURE);

      // Set simulation state
      spVis->setInitialState(spIn);
   }

   void BoussinesqFPlane3DQGModel::addAsciiOutputFiles(SharedSimulation spSim)
   {
      // Create Nusselt number writer
      IoVariable::SharedCartesian1DNusseltZWriter spNusselt(new IoVariable::Cartesian1DNusseltZWriter(SchemeType::type()));
      spNusselt->expect(PhysicalNames::DZ_MEANTEMPERATURE);
      spSim->addAsciiOutputFile(spNusselt);

      // Create temperature energy writer
      IoVariable::SharedCartesian1DScalarEnergyWriter spTemp(new IoVariable::Cartesian1DScalarEnergyWriter("temperature", SchemeType::type()));
      spTemp->expect(PhysicalNames::TEMPERATURE);
      spSim->addAsciiOutputFile(spTemp);

      // Create kinetic energy writer
      IoVariable::SharedCartesian1DStreamEnergyWriter spStream(new IoVariable::Cartesian1DStreamEnergyWriter("kinetic", SchemeType::type()));
      spStream->expect(PhysicalNames::STREAMFUNCTION);
      spStream->expect(PhysicalNames::VELOCITYZ);
      spSim->addAsciiOutputFile(spStream);
   }

   void BoussinesqFPlane3DQGModel::addHdf5OutputFiles(SharedSimulation spSim)
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

      // Add mean temperature to ouput file
      spState->expect(PhysicalNames::DZ_MEANTEMPERATURE);

      spSim->addHdf5OutputFile(spState);
   }

   void BoussinesqFPlane3DQGModel::setInitialState(SharedSimulation spSim)
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
