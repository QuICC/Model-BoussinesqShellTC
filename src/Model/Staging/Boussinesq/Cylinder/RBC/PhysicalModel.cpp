/** 
 * @file PhysicalModel.cpp
 * @brief Source of Boussinesq Rayleigh-Benard convection in a cylinder (toroidal-poloidal formulation)
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
#include MAKE_STR( QUICC_MODEL_PATH/Boussinesq/Cylinder/RBC/PhysicalModel.hpp )

// Project includes
//
#include MAKE_STR( QUICC_MODEL_PATH/Boussinesq/Cylinder/RBC/Transport.hpp )
#include MAKE_STR( QUICC_MODEL_PATH/Boussinesq/Cylinder/RBC/Momentum.hpp )
#include "Enums/FieldIds.hpp"
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
#include "Generator/Visualizers/NonlinearVectorFieldVisualizer.hpp"
#include "Model/PhysicalModelBase.hpp"

namespace QuICC {

namespace Model {

namespace Boussinesq {

namespace Cylinder {

namespace RBC {

   const std::string PhysicalModel::PYMODULE = "boussinesq_rbccylinder";

   const std::string PhysicalModel::PYCLASS = "BoussinesqRBCCylinder";

   void PhysicalModel::addEquations(SharedSimulation spSim)
   {
      // Add transport equation
      spSim->addScalarEquation<Equations::Boussinesq::Cylinder::RBC::Transport>();
                                                                
      // Add Navier-Stokes equation                             
      spSim->addVectorEquation<Equations::Boussinesq::Cylinder::RBC::Momentum>();
   }

   void PhysicalModel::addStates(SharedStateGenerator spGen)
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
         spVector->setSpectralType(FieldComponents::Spectral::TOR, Equations::CylinderExactStateIds::SPEC_UNIT);
         spVector->setSpectralType(FieldComponents::Spectral::POL, Equations::CylinderExactStateIds::SPEC_UNIT);

         // Add scalar exact initial state generator
         spScalar = spGen->addScalarEquation<Equations::CylinderExactScalarState>();
         spScalar->setIdentity(PhysicalNames::TEMPERATURE);
         spScalar->setSpectralType(Equations::CylinderExactStateIds::SPEC_UNIT);

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

   void PhysicalModel::addVisualizers(SharedVisualizationGenerator spVis)
   {
      // Shared pointer to basic field visualizer
      Equations::SharedScalarFieldVisualizer spScalar;
      Equations::SharedScalarFieldTrivialVisualizer spSTrivial;
      Equations::SharedVectorFieldVisualizer spVector;
      Equations::SharedVectorFieldTrivialVisualizer spVTrivial;
      Equations::SharedNonlinearVectorFieldVisualizer spVNL;

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

      // Add output file
      IoVariable::SharedVisualizationFileWriter spOut(new IoVariable::VisualizationFileWriter(SchemeType::type()));
      spOut->expect(PhysicalNames::TEMPERATURE);
      spOut->expect(PhysicalNames::MEAN_TEMPERATURE);
      spOut->expect(PhysicalNames::FLUCT_TEMPERATURE);
      spOut->expect(PhysicalNames::VELOCITY);
      spVis->addHdf5OutputFile(spOut);
   }

   void PhysicalModel::setVisualizationState(SharedVisualizationGenerator spVis)
   {
      // Create and add initial state file to IO
      IoVariable::SharedStateFileReader spIn(new IoVariable::StateFileReader("4Visu", SchemeType::type(), SchemeType::isRegular()));

      // Set expected fields
      spIn->expect(PhysicalNames::TEMPERATURE);
      spIn->expect(PhysicalNames::VELOCITY);

      // Set simulation state
      spVis->setInitialState(spIn);
   }

   void PhysicalModel::addAsciiOutputFiles(SharedSimulation spSim)
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
