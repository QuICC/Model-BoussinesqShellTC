/** 
 * @file Model.cpp
 * @brief Source of the rotating Boussinesq thermal dynamo in a plane layer (toroidal/poloidal formulation) model
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

/// Define small macros allowing to convert to string
#define MAKE_STR_X( _P ) # _P
#define MAKE_STR( _P ) MAKE_STR_X( _P )

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include MAKE_STR( QUICC_MODEL_PATH/Boussinesq/Plane/Dynamo/PhysicalModel.hpp )

// Project includes
//
#include MAKE_STR( QUICC_MODEL_PATH/Boussinesq/Plane/Dynamo/Transport.hpp )
#include MAKE_STR( QUICC_MODEL_PATH/Boussinesq/Plane/Dynamo/Momentum.hpp )
#include MAKE_STR( QUICC_MODEL_PATH/Boussinesq/Plane/Dynamo/Induction.hpp )
#include "Enums/FieldIds.hpp"
#include "IoVariable/StateFileReader.hpp"
#include "IoVariable/StateFileWriter.hpp"
#include "IoVariable/Cartesian1DScalarEnergyWriter.hpp"
#include "IoVariable/Cartesian1DTorPolEnergyWriter.hpp"
#include "IoVariable/Cartesian1DNusseltDZWriter.hpp"
#include "IoVariable/VisualizationFileWriter.hpp"
#include "IoTools/IdToHuman.hpp"
#include "Generator/States/RandomScalarState.hpp"
#include "Generator/States/RandomVectorState.hpp"
#include "Generator/States/CartesianExactScalarState.hpp"
#include "Generator/States/CartesianExactVectorState.hpp"
#include "Generator/Visualizers/ScalarFieldVisualizer.hpp"
#include "Generator/Visualizers/ScalarFieldTrivialVisualizer.hpp"
#include "Generator/Visualizers/VectorFieldVisualizer.hpp"
#include "Model/PhysicalModelBase.hpp"

namespace QuICC {

namespace Model {

namespace Boussinesq {

namespace Plane {

namespace Dynamo {

   const std::string PhysicalModel::PYMODULE = "boussinesq_dynamoplane";

   const std::string PhysicalModel::PYCLASS = "BoussinesqDynamoPlane";

   void PhysicalModel::addEquations(SharedSimulation spSim)
   {
      // Add transport equation
      spSim->addScalarEquation<Equations::Boussinesq::Plane::Dynamo::Transport>();
      
      // Add Navier-Stokes equation
      spSim->addVectorEquation<Equations::Boussinesq::Plane::Dynamo::Momentum>();
      
      // Add Induction equation
      spSim->addVectorEquation<Equations::Boussinesq::Plane::Dynamo::Induction>();
   }

   void PhysicalModel::addStates(SharedStateGenerator spGen)
   {
      // Generate "exact" solutions (trigonometric or monomial)
      if(false)
      {
         // Shared pointer to equation
         Equations::SharedCartesianExactScalarState spScalar;
         Equations::SharedCartesianExactVectorState spVector;
//         Equations::SharedCartesianExactVectorState spExact;

         // Add vector exact initial state generator
         spVector = spGen->addVectorEquation<Equations::CartesianExactVectorState>();
         spVector->setIdentity(PhysicalNames::VELOCITY);
         spVector->setStateType(Equations::CartesianExactStateIds::TORPOLTFF);

         // Add vector exact initial state generator
         spVector = spGen->addVectorEquation<Equations::CartesianExactVectorState>();
         spVector->setIdentity(PhysicalNames::MAGNETIC);
         spVector->setStateType(Equations::CartesianExactStateIds::TORPOLTFF);

         // Add scalar exact initial state generator
         spScalar = spGen->addScalarEquation<Equations::CartesianExactScalarState>();
         spScalar->setIdentity(PhysicalNames::TEMPERATURE);
         spScalar->setStateType(Equations::CartesianExactStateIds::PLANFORMSQUARES);
         spScalar->setModeOptions(1e0, 10.0, 1e0, 10.0);

         // Add imposed magnetic field
//         spExact = spGen->addVectorEquation<Equations::CartesianExactVectorState>();
//         spExact->setIdentity(PhysicalNames::IMPOSED_MAGNETIC);
//         spExact->setStateType(Equations::CartesianExactStateIds::TORPOLCNST);

      // Generate random spectrum
      } else
      {
         // Shared pointer to random initial state equation
         Equations::SharedRandomScalarState spScalar;
         Equations::SharedRandomVectorState spVector;
//         Equations::SharedCartesianExactVectorState spExact;

         // Add scalar random initial state generator 
         spVector = spGen->addVectorEquation<Equations::RandomVectorState>();
         spVector->setIdentity(PhysicalNames::VELOCITY);
         spVector->setSpectrum(FieldComponents::Spectral::TOR, -1e-6, 1e-6, 1e4, 1e4, 1e4);
         spVector->setSpectrum(FieldComponents::Spectral::POL, -1e-6, 1e-6, 1e4, 1e4, 1e4);

         // Add scalar random initial state generator 
         spVector = spGen->addVectorEquation<Equations::RandomVectorState>();
         spVector->setIdentity(PhysicalNames::MAGNETIC);
         spVector->setSpectrum(FieldComponents::Spectral::TOR, -1e-10, 1e-10, 1e4, 1e4, 1e4);
         spVector->setSpectrum(FieldComponents::Spectral::POL, -1e-10, 1e-10, 1e4, 1e4, 1e4);

         // Add imposed magnetic field
//         spExact = spGen->addVectorEquation<Equations::CartesianExactVectorState>();
//         spExact->setIdentity(PhysicalNames::IMPOSED_MAGNETIC);
//         spExact->setStateType(Equations::CartesianExactStateIds::TORPOLCNST);

         // Add scalar random initial state generator
         spScalar = spGen->addScalarEquation<Equations::RandomScalarState>();
         spScalar->setIdentity(PhysicalNames::TEMPERATURE);
         spScalar->setSpectrum(-1e-6, 1e-6, 1e4, 1e4, 1e4);
      }

      // Add output file
      IoVariable::SharedStateFileWriter spOut(new IoVariable::StateFileWriter(SchemeType::type(), SchemeType::isRegular()));
      spOut->expect(PhysicalNames::VELOCITY);
      spOut->expect(PhysicalNames::MAGNETIC);
      spOut->expect(PhysicalNames::TEMPERATURE);
//      spOut->expect(PhysicalNames::IMPOSED_MAGNETIC);
      spGen->addHdf5OutputFile(spOut);
   }

   void PhysicalModel::addVisualizers(SharedVisualizationGenerator spVis)
   {
      // Shared pointer to basic field visualizer
      Equations::SharedScalarFieldVisualizer spScalar;
      Equations::SharedScalarFieldTrivialVisualizer spSTrivial;
      Equations::SharedVectorFieldVisualizer spVector;

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

      // Add velocity field visualization
      spVector = spVis->addVectorEquation<Equations::VectorFieldVisualizer>();
      spVector->setFields(true, false, true);
      spVector->setIdentity(PhysicalNames::VELOCITY);

      // Add magnetic field visualization
      spVector = spVis->addVectorEquation<Equations::VectorFieldVisualizer>();
      spVector->setFields(true, false, false);
      spVector->setIdentity(PhysicalNames::MAGNETIC);

      // Add imposed magnetic field visualization
//      spVector = spVis->addVectorEquation<Equations::VectorFieldVisualizer>();
//      spVector->setFields(true, false, false);
//      spVector->setIdentity(PhysicalNames::IMPOSED_MAGNETIC);


      // Add output file
      IoVariable::SharedVisualizationFileWriter spOut(new IoVariable::VisualizationFileWriter(SchemeType::type()));
      spOut->expect(PhysicalNames::TEMPERATURE);
      spOut->expect(PhysicalNames::MEAN_TEMPERATURE);
      spOut->expect(PhysicalNames::FLUCT_TEMPERATURE);
      spOut->expect(PhysicalNames::VELOCITY);
      spOut->expect(PhysicalNames::MAGNETIC);
//      spOut->expect(PhysicalNames::IMPOSED_MAGNETIC);
      spVis->addHdf5OutputFile(spOut);
   }

   void PhysicalModel::setVisualizationState(SharedVisualizationGenerator spVis)
   {
      // Create and add initial state file to IO
      IoVariable::SharedStateFileReader spIn(new IoVariable::StateFileReader("4Visu", SchemeType::type(), SchemeType::isRegular()));

      // Set expected fields
      spIn->expect(PhysicalNames::TEMPERATURE);
      spIn->expect(PhysicalNames::VELOCITY);
      spIn->expect(PhysicalNames::MAGNETIC);
//      spIn->expect(PhysicalNames::IMPOSED_MAGNETIC);
      // Set simulation state
      spVis->setInitialState(spIn);
   }

   void PhysicalModel::addAsciiOutputFiles(SharedSimulation spSim)
   {
      // Create temperature energy writer
      IoVariable::SharedCartesian1DScalarEnergyWriter spTemp(new IoVariable::Cartesian1DScalarEnergyWriter("temperature", SchemeType::type()));
      spTemp->expect(PhysicalNames::TEMPERATURE);
      spSim->addAsciiOutputFile(spTemp);

      // Create kinetic energy writer
      IoVariable::SharedCartesian1DTorPolEnergyWriter spKinetic(new IoVariable::Cartesian1DTorPolEnergyWriter("kinetic", SchemeType::type()));
      spKinetic->expect(PhysicalNames::VELOCITY);
      spSim->addAsciiOutputFile(spKinetic);

      // Create magnetic energy writer
      IoVariable::SharedCartesian1DTorPolEnergyWriter spMagnetic(new IoVariable::Cartesian1DTorPolEnergyWriter("magnetic", SchemeType::type()));
      spMagnetic->expect(PhysicalNames::MAGNETIC);
      spSim->addAsciiOutputFile(spMagnetic);

      // Create nusselt number writer
      IoVariable::SharedCartesian1DNusseltDZWriter spNusselt(new IoVariable::Cartesian1DNusseltDZWriter(SchemeType::type()));
      spNusselt->expect(PhysicalNames::TEMPERATURE);
      spSim->addAsciiOutputFile(spNusselt);
   }

   void PhysicalModel::addHdf5OutputFiles(SharedSimulation spSim)
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

   void PhysicalModel::addStatsOutputFiles(SharedSimulation spSim)
   {
   }

   void PhysicalModel::setInitialState(SharedSimulation spSim)
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
}
}
}
}
