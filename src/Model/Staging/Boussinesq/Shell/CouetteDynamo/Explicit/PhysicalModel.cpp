/** 
 * @file PhysicalModel.cpp
 * @brief Source of the Boussinesq spherical Couette dynamo in a spherical shell (Toroidal/Poloidal formulation) without coupled solve (standard implementation)
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
#include MAKE_STR( QUICC_MODEL_PATH/Boussinesq/Shell/CouetteDynamo/Explicit/PhysicalModel.hpp )

// Project includes
//
#include MAKE_STR( QUICC_MODEL_PATH/Boussinesq/Shell/Couette/Momentum.hpp )
#include MAKE_STR( QUICC_MODEL_PATH/Boussinesq/Shell/CouetteDynamo/Momentum.hpp )
#include MAKE_STR( QUICC_MODEL_PATH/Boussinesq/Shell/CouetteDynamo/Induction.hpp )
#include "Enums/FieldIds.hpp"
#include "IoVariable/StateFileReader.hpp"
#include "IoVariable/StateFileWriter.hpp"
#include "IoVariable/VisualizationFileWriter.hpp"
#include "IoTools/IdToHuman.hpp"
#include "IoVariable/ShellTorPolEnergyWriter.hpp"
#include "Generator/States/RandomVectorState.hpp"
#include "Generator/States/ShellExactStateIds.hpp"
#include "Generator/States/ShellExactVectorState.hpp"
#include "Generator/Visualizers/VectorFieldVisualizer.hpp"
#include "Generator/Visualizers/SphericalVerticalFieldVisualizer.hpp"
#include "Model/PhysicalModelBase.hpp"

namespace QuICC {

namespace Model {

namespace Boussinesq {

namespace Shell {

namespace CouetteDynamo {

namespace Explicit {

   const std::string PhysicalModel::PYMODULE = "boussinesq_dynamocouetteshell_std";

   const std::string PhysicalModel::PYCLASS = "BoussinesqDynamoCouetteShellStd";

   void PhysicalModel::addEquations(SharedSimulation spSim)
   {
      // Add Navier-Stokes equation for kinematic dynamo
      spSim->addVectorEquation<Equations::Boussinesq::Shell::Couette::Momentum>();
      // Add Navier-Stokes equation for fully nonlinear dynamo
      //spSim->addVectorEquation<Equations::Boussinesq::Shell::CouetteDynamo::Momentum>();

      // Add induction equation
      spSim->addVectorEquation<Equations::Boussinesq::Shell::CouetteDynamo::Induction>();
   }

   void PhysicalModel::addStates(SharedStateGenerator spGen)
   {
      // Generate "exact" solutions (trigonometric or monomial)
      if(true)
      {
         // Shared pointer to equation
         Equations::SharedShellExactVectorState spVector;

         std::vector<std::tr1::tuple<int,int,MHDComplex> > tSH;

         // Add velocity initial state generator
         spVector = spGen->addVectorEquation<Equations::ShellExactVectorState>();
         spVector->setIdentity(PhysicalNames::VELOCITY);
         switch(4)
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

            case 3:
               spVector->setStateType(Equations::ShellExactStateIds::BENCHVELC1);
               break;

            case 4:
               spVector->setStateType(Equations::ShellExactStateIds::NOISE);
               break;
         }

         // Add magnetic field initial state generator
         spVector = spGen->addVectorEquation<Equations::ShellExactVectorState>();
         spVector->setIdentity(PhysicalNames::MAGNETIC);
         switch(4)
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

            case 3:
               spVector->setStateType(Equations::ShellExactStateIds::BENCHMAGC1);
               break;

            case 4:
               spVector->setStateType(Equations::ShellExactStateIds::NOISE);
               break;
         }

         // Add imposed magnetic field initial state generator
         spVector = spGen->addVectorEquation<Equations::ShellExactVectorState>();
         spVector->setIdentity(PhysicalNames::IMPOSED_MAGNETIC);
         switch(3)
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

            case 3:
               spVector->setStateType(Equations::ShellExactStateIds::IMPMAGMATSUI);
               break;

            case 4:
               spVector->setStateType(Equations::ShellExactStateIds::NOISE);
               break;
         }

      // Generate random spectrum
      } else
      {
         // Shared pointer to random initial state equation
         Equations::SharedRandomVectorState spVector;

         // Add vector random initial state generator 
         spVector = spGen->addVectorEquation<Equations::RandomVectorState>();
         spVector->setIdentity(PhysicalNames::VELOCITY);
         spVector->setSpectrum(FieldComponents::Spectral::TOR, -1e-4, 1e-4, 1e4, 1e4, 1e4);
         spVector->setSpectrum(FieldComponents::Spectral::POL, -1e-4, 1e-4, 1e4, 1e4, 1e4);

         // Add vector random initial state generator 
         spVector = spGen->addVectorEquation<Equations::RandomVectorState>();
         spVector->setIdentity(PhysicalNames::MAGNETIC);
         spVector->setSpectrum(FieldComponents::Spectral::TOR, -1e-4, 1e-4, 1e4, 1e4, 1e4);
         spVector->setSpectrum(FieldComponents::Spectral::POL, -1e-4, 1e-4, 1e4, 1e4, 1e4);

         // Add vector random initial state generator 
         spVector = spGen->addVectorEquation<Equations::RandomVectorState>();
         spVector->setIdentity(PhysicalNames::IMPOSED_MAGNETIC);
         spVector->setSpectrum(FieldComponents::Spectral::TOR, -1e-4, 1e-4, 1e4, 1e4, 1e4);
         spVector->setSpectrum(FieldComponents::Spectral::POL, -1e-4, 1e-4, 1e4, 1e4, 1e4);
      }

      // Add output file
      IoVariable::SharedStateFileWriter spOut(new IoVariable::StateFileWriter(SchemeType::type(), SchemeType::isRegular()));
      spOut->expect(PhysicalNames::VELOCITY);
      spOut->expect(PhysicalNames::MAGNETIC);
      spOut->expect(PhysicalNames::IMPOSED_MAGNETIC);
      spGen->addHdf5OutputFile(spOut);
   }

   void PhysicalModel::addVisualizers(SharedVisualizationGenerator spVis)
   {
      // Shared pointer to basic field visualizer
      Equations::SharedVectorFieldVisualizer spVector;
      Equations::SharedSphericalVerticalFieldVisualizer spVertical;

      // Add velocity field visualization
      spVector = spVis->addVectorEquation<Equations::VectorFieldVisualizer>();
      spVector->setFields(true, false, false);
      spVector->setIdentity(PhysicalNames::VELOCITY);

      // Add vertical velocity visualization
      spVertical = spVis->addScalarEquation<Equations::SphericalVerticalFieldVisualizer>();
      spVertical->setFieldType(FieldType::VECTOR);
      spVertical->setIdentity(PhysicalNames::VELOCITYZ, PhysicalNames::VELOCITY);

      // Add vertical vorticity visualization
      spVertical = spVis->addScalarEquation<Equations::SphericalVerticalFieldVisualizer>();
      spVertical->setFieldType(FieldType::CURL);
      spVertical->setIdentity(PhysicalNames::VORTICITYZ, PhysicalNames::VELOCITY);

      // Add magnetic field visualization
      spVector = spVis->addVectorEquation<Equations::VectorFieldVisualizer>();
      spVector->setFields(true, false, false);
      spVector->setIdentity(PhysicalNames::MAGNETIC);

      // Add imposed magnetic field visualization
      spVector = spVis->addVectorEquation<Equations::VectorFieldVisualizer>();
      spVector->setFields(true, false, false);
      spVector->setIdentity(PhysicalNames::IMPOSED_MAGNETIC);

      // Add output file
      IoVariable::SharedVisualizationFileWriter spOut(new IoVariable::VisualizationFileWriter(SchemeType::type()));
      spOut->expect(PhysicalNames::VELOCITY);
      spOut->expect(PhysicalNames::VELOCITYZ);
      spOut->expect(PhysicalNames::VORTICITYZ);
      spOut->expect(PhysicalNames::MAGNETIC);
      spOut->expect(PhysicalNames::IMPOSED_MAGNETIC);
      spVis->addHdf5OutputFile(spOut);
   }

   void PhysicalModel::setVisualizationState(SharedVisualizationGenerator spVis)
   {
      // Create and add initial state file to IO
      IoVariable::SharedStateFileReader spIn(new IoVariable::StateFileReader("4Visu", SchemeType::type(), SchemeType::isRegular()));

      // Set expected fields
      spIn->expect(PhysicalNames::VELOCITY);
      spIn->expect(PhysicalNames::MAGNETIC);
      spIn->expect(PhysicalNames::IMPOSED_MAGNETIC);

      // Set simulation state
      spVis->setInitialState(spIn);
   }

   void PhysicalModel::addAsciiOutputFiles(SharedSimulation spSim)
   {
      // Create kinetic energy writer
      IoVariable::SharedShellTorPolEnergyWriter spVector(new IoVariable::ShellTorPolEnergyWriter("kinetic", SchemeType::type()));
      spVector->expect(PhysicalNames::VELOCITY);
      spSim->addAsciiOutputFile(spVector);

      // Create magnetic energy writer
      spVector = IoVariable::SharedShellTorPolEnergyWriter(new IoVariable::ShellTorPolEnergyWriter("magnetic", SchemeType::type()));
      spVector->expect(PhysicalNames::MAGNETIC);
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
      spState->expect(PhysicalNames::IMPOSED_MAGNETIC);
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
      spInit->expect(PhysicalNames::IMPOSED_MAGNETIC);

      // Set simulation state
      spSim->setInitialState(spInit);
   }

}
}
}
}
}
}
