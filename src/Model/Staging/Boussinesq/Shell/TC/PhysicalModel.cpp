/** 
 * @file PhysicalModel.cpp
 * @brief Source of the Boussinesq thermal convection in a spherical shell (Toroidal/Poloidal formulation)
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
#include MAKE_STR( QUICC_MODEL_PATH/Boussinesq/Shell/TC/PhysicalModel.hpp )

// Project includes
//
#include MAKE_STR( QUICC_MODEL_PATH/Boussinesq/Shell/TC/Transport.hpp )
#include MAKE_STR( QUICC_MODEL_PATH/Boussinesq/Shell/TC/Momentum.hpp )
#include "Enums/FieldIds.hpp"
#include "IoVariable/StateFileReader.hpp"
#include "IoVariable/StateFileWriter.hpp"
#include "IoVariable/VisualizationFileWriter.hpp"
#include "IoTools/IdToHuman.hpp"
#include "IoVariable/ShellScalarEnergyWriter.hpp"
#include "IoVariable/ShellTorPolEnergyWriter.hpp"
#include "Generator/States/RandomScalarState.hpp"
#include "Generator/States/RandomVectorState.hpp"
#include "Generator/States/ShellExactStateIds.hpp"
#include "Generator/States/ShellExactScalarState.hpp"
#include "Generator/States/ShellExactVectorState.hpp"
#include "Generator/Visualizers/ScalarFieldVisualizer.hpp"
#include "Generator/Visualizers/VectorFieldVisualizer.hpp"
#include "Model/PhysicalModelBase.hpp"

namespace QuICC {

namespace Model {

namespace Boussinesq {

namespace Shell {

namespace TC {

   const std::string PhysicalModel::PYMODULE = "boussinesq_tcshell";

   const std::string PhysicalModel::PYCLASS = "BoussinesqTCShell";

   void PhysicalModel::addEquations(SharedSimulation spSim)
   {
      // Add transport equation
      spSim->addScalarEquation<Equations::Boussinesq::Shell::TC::Transport>();
                                                           
      // Add Navier-Stokes equation                        
      spSim->addVectorEquation<Equations::Boussinesq::Shell::TC::Momentum>();
   }

   void PhysicalModel::addStates(SharedStateGenerator spGen)
   {
      // Generate "exact" solutions (trigonometric or monomial)
      if(true)
      {
         // Shared pointer to equation
         Equations::SharedShellExactScalarState spScalar;
         Equations::SharedShellExactVectorState spVector;

         std::vector<std::tr1::tuple<int,int,MHDComplex> > tSH;

         // Add temperature initial state generator
         spScalar = spGen->addScalarEquation<Equations::ShellExactScalarState>();
         spScalar->setIdentity(PhysicalNames::TEMPERATURE);
         spScalar->setStateType(Equations::ShellExactStateIds::HARMONIC);
         tSH.clear(); 
         tSH.push_back(std::tr1::make_tuple(0,0,MHDComplex(1,1)));
         tSH.push_back(std::tr1::make_tuple(1,0,MHDComplex(1,1)));
         tSH.push_back(std::tr1::make_tuple(1,1,MHDComplex(1,1)));
         tSH.push_back(std::tr1::make_tuple(2,0,MHDComplex(1,1)));
         tSH.push_back(std::tr1::make_tuple(2,1,MHDComplex(1,1)));
         tSH.push_back(std::tr1::make_tuple(2,2,MHDComplex(1,1)));
         tSH.push_back(std::tr1::make_tuple(5,5,MHDComplex(1,1)));
         spScalar->setHarmonicOptions(tSH);
         spScalar->setHarmonicOptions(tSH);

         // Add velocity initial state generator
         spVector = spGen->addVectorEquation<Equations::ShellExactVectorState>();
         spVector->setIdentity(PhysicalNames::VELOCITY);
         switch(2)
         {
            case 0:
               spVector->setStateType(Equations::ShellExactStateIds::TOROIDAL);
               tSH.clear(); 
               tSH.push_back(std::tr1::make_tuple(0,0,MHDComplex(1,0)));
               tSH.push_back(std::tr1::make_tuple(1,0,MHDComplex(1,0)));
               tSH.push_back(std::tr1::make_tuple(1,1,MHDComplex(1,0)));
               tSH.push_back(std::tr1::make_tuple(2,0,MHDComplex(1,0)));
               tSH.push_back(std::tr1::make_tuple(2,1,MHDComplex(1,0)));
               tSH.push_back(std::tr1::make_tuple(2,2,MHDComplex(1,0)));
               tSH.push_back(std::tr1::make_tuple(5,4,MHDComplex(1,0)));
               spVector->setHarmonicOptions(FieldComponents::Spectral::TOR, tSH);
               break;

            case 1:
               spVector->setStateType(Equations::ShellExactStateIds::POLOIDAL);
               tSH.clear(); 
               tSH.push_back(std::tr1::make_tuple(0,0,MHDComplex(1,0)));
               tSH.push_back(std::tr1::make_tuple(1,0,MHDComplex(1,0)));
               tSH.push_back(std::tr1::make_tuple(1,1,MHDComplex(1,0)));
               tSH.push_back(std::tr1::make_tuple(2,0,MHDComplex(1,0)));
               tSH.push_back(std::tr1::make_tuple(2,1,MHDComplex(1,0)));
               tSH.push_back(std::tr1::make_tuple(2,2,MHDComplex(1,0)));
               tSH.push_back(std::tr1::make_tuple(4,3,MHDComplex(1,0)));
               spVector->setHarmonicOptions(FieldComponents::Spectral::POL, tSH);
               break;

            case 2:
               spVector->setStateType(Equations::ShellExactStateIds::TORPOL);
               tSH.clear(); 
               tSH.push_back(std::tr1::make_tuple(0,0,MHDComplex(1,0)));
               tSH.push_back(std::tr1::make_tuple(1,0,MHDComplex(1,0)));
               tSH.push_back(std::tr1::make_tuple(1,1,MHDComplex(1,0)));
               tSH.push_back(std::tr1::make_tuple(2,0,MHDComplex(1,0)));
               tSH.push_back(std::tr1::make_tuple(2,1,MHDComplex(1,0)));
               tSH.push_back(std::tr1::make_tuple(2,2,MHDComplex(1,0)));
               tSH.push_back(std::tr1::make_tuple(5,4,MHDComplex(1,0)));
               spVector->setHarmonicOptions(FieldComponents::Spectral::TOR, tSH);
               tSH.clear(); 
               tSH.push_back(std::tr1::make_tuple(0,0,MHDComplex(1,0)));
               tSH.push_back(std::tr1::make_tuple(1,0,MHDComplex(1,0)));
               tSH.push_back(std::tr1::make_tuple(1,1,MHDComplex(1,0)));
               tSH.push_back(std::tr1::make_tuple(2,0,MHDComplex(1,0)));
               tSH.push_back(std::tr1::make_tuple(2,1,MHDComplex(1,0)));
               tSH.push_back(std::tr1::make_tuple(2,2,MHDComplex(1,0)));
               tSH.push_back(std::tr1::make_tuple(4,3,MHDComplex(1,0)));
               spVector->setHarmonicOptions(FieldComponents::Spectral::POL, tSH);
               break;
         }

      // Generate random spectrum
      } else
      {
         // Shared pointer to random initial state equation
         Equations::SharedRandomScalarState spScalar;
         Equations::SharedRandomVectorState spVector;

         // Add scalar random initial state generator 
         spVector = spGen->addVectorEquation<Equations::RandomVectorState>();
         spVector->setIdentity(PhysicalNames::VELOCITY);
         spVector->setSpectrum(FieldComponents::Spectral::TOR, -1e-2, 1e-2, 1e4, 1e4, 1e4);
         spVector->setSpectrum(FieldComponents::Spectral::POL, -1e-2, 1e-2, 1e4, 1e4, 1e4);

         // Add scalar random initial state generator
         spScalar = spGen->addScalarEquation<Equations::RandomScalarState>();
         spScalar->setIdentity(PhysicalNames::TEMPERATURE);
         spScalar->setSpectrum(-1e-1, 1e-1, 1e4, 1e4, 1e4);
      }

      // Add output file
      IoVariable::SharedStateFileWriter spOut(new IoVariable::StateFileWriter(SchemeType::type(), SchemeType::isRegular()));
      spOut->expect(PhysicalNames::TEMPERATURE);
      spOut->expect(PhysicalNames::VELOCITY);
      spGen->addHdf5OutputFile(spOut);
   }

   void PhysicalModel::addVisualizers(SharedVisualizationGenerator spVis)
   {
      // Shared pointer to basic field visualizer
      Equations::SharedScalarFieldVisualizer spScalar;
      Equations::SharedVectorFieldVisualizer spVector;

      // Add temperature field visualization
      spScalar = spVis->addScalarEquation<Equations::ScalarFieldVisualizer>();
      spScalar->setFields(true, true);
      spScalar->setIdentity(PhysicalNames::TEMPERATURE);

      // Add velocity field visualization
      spVector = spVis->addVectorEquation<Equations::VectorFieldVisualizer>();
      spVector->setFields(true, false, true);
      spVector->setIdentity(PhysicalNames::VELOCITY);

      // Add output file
      IoVariable::SharedVisualizationFileWriter spOut(new IoVariable::VisualizationFileWriter(SchemeType::type()));
      spOut->expect(PhysicalNames::TEMPERATURE);
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
      IoVariable::SharedShellScalarEnergyWriter spScalar(new IoVariable::ShellScalarEnergyWriter("temperature", SchemeType::type()));
      spScalar->expect(PhysicalNames::TEMPERATURE);
      spSim->addAsciiOutputFile(spScalar);

      // Create kinetic energy writer
      IoVariable::SharedShellTorPolEnergyWriter spVector(new IoVariable::ShellTorPolEnergyWriter("kinetic", SchemeType::type()));
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
