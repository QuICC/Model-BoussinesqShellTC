/** 
 * @file PhysicalModel.cpp
 * @brief Source of the Boussinesq precession in a sphere (Toroidal/Poloidal formulation) without coupled solve (standard implementation)
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
#include MAKE_STR( QUICC_MODEL_PATH/Boussinesq/Sphere/Precession/Explicit/PhysicalModel.hpp )

// Project includes
//
#include MAKE_STR( QUICC_MODEL_PATH/Boussinesq/Sphere/Precession/Momentum.hpp )
#include "Enums/FieldIds.hpp"
#include "IoVariable/StateFileReader.hpp"
#include "IoVariable/StateFileWriter.hpp"
#include "IoVariable/VisualizationFileWriter.hpp"
#include "IoTools/IdToHuman.hpp"
#include "IoVariable/SphereTorPolEnergyWriter.hpp"
#include "Generator/States/RandomVectorState.hpp"
#include "Generator/States/SphereExactStateIds.hpp"
#include "Generator/States/SphereExactVectorState.hpp"
#include "Generator/Visualizers/VectorFieldVisualizer.hpp"
#include "Generator/Visualizers/SphericalVerticalFieldVisualizer.hpp"
#include "Model/PhysicalModelBase.hpp"

namespace QuICC {

namespace Model {

namespace Boussinesq {

namespace Sphere {

namespace Precession {

namespace Explicit {

   const std::string PhysicalModel::PYMODULE = "boussinesq_precessionsphere_std";

   const std::string PhysicalModel::PYCLASS = "BoussinesqPrecessionSphereStd";

   void PhysicalModel::addEquations(SharedSimulation spSim)
   {
      // Add Navier-Stokes equation
      spSim->addVectorEquation<Equations::Boussinesq::Sphere::Precession::Momentum>();
   }

   void PhysicalModel::addStates(SharedStateGenerator spGen)
   {
      // Generate "exact" solutions (trigonometric or monomial)
      if(true)
      {
         // Shared pointer to equation
         Equations::SharedSphereExactVectorState spVector;

         Equations::SphereExactVectorState::HarmonicModeType tSH;
         std::pair<Equations::SphereExactVectorState::HarmonicModeType::iterator,bool> ptSH; 

         // Add velocity initial state generator
         spVector = spGen->addVectorEquation<Equations::SphereExactVectorState>();
         spVector->setIdentity(PhysicalNames::VELOCITY);
         switch(3)
         {
            case 0:
               // Toroidal
               spVector->setSpectralType(Equations::SphereExactStateIds::HARMONIC);
               tSH.clear(); 
               ptSH = tSH.insert(std::make_pair(std::make_pair(1,1), std::map<int,MHDComplex>()));
               ptSH.first->second.insert(std::make_pair(7, MHDComplex(1.0)));
               spVector->setHarmonicOptions(FieldComponents::Spectral::TOR, tSH);
               break;

            case 1:
               // Poloidal
               spVector->setSpectralType(Equations::SphereExactStateIds::HARMONIC);
               tSH.clear(); 
               ptSH = tSH.insert(std::make_pair(std::make_pair(2,0), std::map<int,MHDComplex>()));
               ptSH.first->second.insert(std::make_pair(7, MHDComplex(1.0)));
               spVector->setHarmonicOptions(FieldComponents::Spectral::POL, tSH);
               break;

            case 2:
               // Toroidal
               spVector->setSpectralType(Equations::SphereExactStateIds::HARMONIC);
               tSH.clear(); 
               ptSH = tSH.insert(std::make_pair(std::make_pair(1,1), std::map<int,MHDComplex>()));
               ptSH.first->second.insert(std::make_pair(7, MHDComplex(1.0)));
               spVector->setHarmonicOptions(FieldComponents::Spectral::TOR, tSH);
               // Poloidal
               tSH.clear(); 
               ptSH = tSH.insert(std::make_pair(std::make_pair(2,0), std::map<int,MHDComplex>()));
               ptSH.first->second.insert(std::make_pair(7, MHDComplex(1.0)));
               spVector->setHarmonicOptions(FieldComponents::Spectral::POL, tSH);
               break;

            case 3:
               spVector->setStateType(Equations::SphereExactStateIds::BENCHVELC1);
               break;

            case 4:
               spVector->setStateType(Equations::SphereExactStateIds::BENCHVELC2);
               break;
         }

      // Generate random spectrum
      } else
      {
         // Shared pointer to random initial state equation
         Equations::SharedRandomVectorState spVector;

         // Add scalar random initial state generator 
         spVector = spGen->addVectorEquation<Equations::RandomVectorState>();
         spVector->setIdentity(PhysicalNames::VELOCITY);
         spVector->setSpectrum(FieldComponents::Spectral::TOR, -1e-4, 1e-4, 1e4, 1e4, 1e4);
         spVector->setSpectrum(FieldComponents::Spectral::POL, -1e-4, 1e-4, 1e4, 1e4, 1e4);
      }

      // Add output file
      IoVariable::SharedStateFileWriter spOut(new IoVariable::StateFileWriter(SchemeType::type(), SchemeType::isRegular()));
      spOut->expect(PhysicalNames::VELOCITY);
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

      // Add output file
      IoVariable::SharedVisualizationFileWriter spOut(new IoVariable::VisualizationFileWriter(SchemeType::type()));
      spOut->expect(PhysicalNames::VELOCITY);
      spOut->expect(PhysicalNames::VELOCITYZ);
      spOut->expect(PhysicalNames::VORTICITYZ);
      spVis->addHdf5OutputFile(spOut);
   }

   void PhysicalModel::setVisualizationState(SharedVisualizationGenerator spVis)
   {
      // Create and add initial state file to IO
      IoVariable::SharedStateFileReader spIn(new IoVariable::StateFileReader("4Visu", SchemeType::type(), SchemeType::isRegular()));

      // Set expected fields
      spIn->expect(PhysicalNames::VELOCITY);

      // Set simulation state
      spVis->setInitialState(spIn);
   }

   void PhysicalModel::addAsciiOutputFiles(SharedSimulation spSim)
   {
      // Create kinetic energy writer
      IoVariable::SharedSphereTorPolEnergyWriter spVector(new IoVariable::SphereTorPolEnergyWriter("kinetic", SchemeType::type()));
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
}
