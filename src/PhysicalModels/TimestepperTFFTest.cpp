/** 
 * @file TimestepperTFFTest.cpp
 * @brief Source of the timestepper test for TFF scheme
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
#include "PhysicalModels/TimestepperTFFTest.hpp"

// Project includes
//
#include "Enums/FieldIds.hpp"
#include "IoVariable/StateFileReader.hpp"
#include "IoVariable/StateFileWriter.hpp"
#include "IoVariable/VisualizationFileWriter.hpp"
#include "IoTools/IdToHuman.hpp"
#include "Equations/Tests/TimestepperTFFScalar.hpp"
#include "Equations/Tests/TimestepperTFFTimeAvgScalar.hpp"
#include "IoVariable/Cartesian1DNusseltZWriter.hpp"
#include "IoVariable/Cartesian1DScalarEnergyWriter.hpp"
#include "IoVariable/Cartesian1DStreamEnergyWriter.hpp"
#include "Generator/States/RandomScalarState.hpp"
#include "Generator/States/CartesianExactScalarState.hpp"
#include "Generator/Visualizers/ScalarFieldVisualizer.hpp"
#include "PhysicalModels/PhysicalModelBase.hpp"

namespace GeoMHDiSCC {

   const std::string TimestepperTFFTest::PYMODULE = "timestepper_tff";

   const std::string TimestepperTFFTest::PYCLASS = "TimestepperTFF";

   void TimestepperTFFTest::addEquations(SharedSimulation spSim)
   {
      // Add timestepper test equation
      spSim->addScalarEquation<Equations::TimestepperTFFScalar>();

      // Add timestepper time averaged test equation
      spSim->addScalarEquation<Equations::TimestepperTFFTimeAvgScalar>();
   }

   void TimestepperTFFTest::addStates(SharedStateGenerator spGen)
   {
      // Generate "exact" solutions (trigonometric or monomial)
      if(true)
      {
         // Shared pointer to equation
         Equations::SharedCartesianExactScalarState spExact;

         // Add initial field
         spExact = spGen->addScalarEquation<Equations::CartesianExactScalarState>();
         spExact->setIdentity(PhysicalNames::TEMPERATURE);
         spExact->setStateType(Equations::CartesianExactStateIds::PEYRET1DA);

         // Add initial field
         spExact = spGen->addScalarEquation<Equations::CartesianExactScalarState>();
         spExact->setIdentity(PhysicalNames::PRESSURE);
         spExact->setStateType(Equations::CartesianExactStateIds::PEYRET1DA);

      // Generate random spectrum
      } else
      {
         // Shared pointer to random initial state equation
         Equations::SharedRandomScalarState spRand;
         Equations::SharedCartesianExactScalarState spExact;

         // Add transport initial state generation equation
         spRand = spGen->addScalarEquation<Equations::RandomScalarState>();
         spRand->setIdentity(PhysicalNames::TEMPERATURE);
         spRand->setSpectrum(-1e-7, 1e-7, 1e4, 1e4, 1e4);

         // Add transport initial state generation equation
         spRand = spGen->addScalarEquation<Equations::RandomScalarState>();
         spRand->setIdentity(PhysicalNames::PRESSURE);
         spRand->setSpectrum(-1e-7, 1e-7, 1e4, 1e4, 1e4);
      }

      // Add output file
      IoVariable::SharedStateFileWriter spOut(new IoVariable::StateFileWriter(SchemeType::type(), SchemeType::isRegular()));
      spOut->expect(PhysicalNames::TEMPERATURE);
      spOut->expect(PhysicalNames::PRESSURE);
      spGen->addHdf5OutputFile(spOut);
   }

   void TimestepperTFFTest::addVisualizers(SharedVisualizationGenerator spVis)
   {
      // Shared pointer to basic field visualizer
      Equations::SharedScalarFieldVisualizer spField;

      // Add transport field visualization
      spField = spVis->addScalarEquation<Equations::ScalarFieldVisualizer>();
      spField->setFields(true, false, false);
      spField->setIdentity(PhysicalNames::TEMPERATURE);

      // Add transport field visualization
      spField = spVis->addScalarEquation<Equations::ScalarFieldVisualizer>();
      spField->setFields(true, false, false);
      spField->setIdentity(PhysicalNames::PRESSURE);

      // Add output file
      IoVariable::SharedVisualizationFileWriter spOut(new IoVariable::VisualizationFileWriter(SchemeType::type()));
      spOut->expect(PhysicalNames::TEMPERATURE);
      spOut->expect(PhysicalNames::PRESSURE);
      spVis->addHdf5OutputFile(spOut);
   }

   void TimestepperTFFTest::setVisualizationState(SharedVisualizationGenerator spVis)
   {
      // Create and add initial state file to IO
      IoVariable::SharedStateFileReader spIn(new IoVariable::StateFileReader("4Visu", SchemeType::type(), SchemeType::isRegular()));

      // Set expected fields
      spIn->expect(PhysicalNames::TEMPERATURE);
      spIn->expect(PhysicalNames::PRESSURE);

      // Set simulation state
      spVis->setInitialState(spIn);
   }

   void TimestepperTFFTest::addAsciiOutputFiles(SharedSimulation spSim)
   {
      // Create temperature energy writer
      IoVariable::SharedCartesian1DScalarEnergyWriter spTemp(new IoVariable::Cartesian1DScalarEnergyWriter("temperature", SchemeType::type()));
      spTemp->expect(PhysicalNames::TEMPERATURE);
      spSim->addAsciiOutputFile(spTemp);

      // Create temperature energy writer
      IoVariable::SharedCartesian1DScalarEnergyWriter spAvgTemp(new IoVariable::Cartesian1DScalarEnergyWriter("time_avg_temperature", SchemeType::type()));
      spAvgTemp->expect(PhysicalNames::PRESSURE);
      spSim->addAsciiOutputFile(spAvgTemp);
   }

   void TimestepperTFFTest::addStatsOutputFiles(SharedSimulation spSim)
   {

   }

   void TimestepperTFFTest::addHdf5OutputFiles(SharedSimulation spSim)
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

   void TimestepperTFFTest::setInitialState(SharedSimulation spSim)
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
      spInit->expect(PhysicalNames::PRESSURE);

      // Set simulation state
      spSim->setInitialState(spInit);
   }

}