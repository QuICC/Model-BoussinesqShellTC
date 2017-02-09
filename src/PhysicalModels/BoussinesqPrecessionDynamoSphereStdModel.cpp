/** 
 * @file BoussinesqPrecessionDynamoSphereStdModel.cpp
 * @brief Source of the Boussinesq precession dynamo in a sphere (Toroidal/Poloidal formulation) without coupled solve (standard implementation)
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
#include "PhysicalModels/BoussinesqPrecessionDynamoSphereStdModel.hpp"

// Project includes
//
#include "Enums/FieldIds.hpp"
#include "IoVariable/StateFileReader.hpp"
#include "IoVariable/StateFileWriter.hpp"
#include "IoVariable/VisualizationFileWriter.hpp"
#include "IoTools/IdToHuman.hpp"
#include "Equations/Sphere/Boussinesq/BoussinesqPrecessionDynamoSphereMomentum.hpp"
#include "Equations/Sphere/Boussinesq/BoussinesqPrecessionDynamoSphereInduction.hpp"
#include "IoVariable/SphereTorPolEnergyWriter.hpp"
#include "Generator/States/RandomVectorState.hpp"
#include "Generator/States/SphereExactStateIds.hpp"
#include "Generator/States/SphereExactVectorState.hpp"
#include "Generator/Visualizers/VectorFieldVisualizer.hpp"
#include "Generator/Visualizers/SphericalVerticalFieldVisualizer.hpp"
#include "PhysicalModels/PhysicalModelBase.hpp"

namespace QuICC {

   const std::string BoussinesqPrecessionDynamoSphereStdModel::PYMODULE = "boussinesq_precessiondynamosphere_std";

   const std::string BoussinesqPrecessionDynamoSphereStdModel::PYCLASS = "BoussinesqPrecessionDynamoSphereStd";

   void BoussinesqPrecessionDynamoSphereStdModel::addEquations(SharedSimulation spSim)
   {
      // Add Navier-Stokes equation
      spSim->addVectorEquation<Equations::BoussinesqPrecessionDynamoSphereMomentum>();
      
      // Add induction equation
      spSim->addVectorEquation<Equations::BoussinesqPrecessionDynamoSphereInduction>();
   }

   void BoussinesqPrecessionDynamoSphereStdModel::addStates(SharedStateGenerator spGen)
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
               tSH.clear(); 
               // Poloidal
               ptSH = tSH.insert(std::make_pair(std::make_pair(2,0), std::map<int,MHDComplex>()));
               ptSH.first->second.insert(std::make_pair(7, MHDComplex(1.0)));
               spVector->setHarmonicOptions(FieldComponents::Spectral::POL, tSH);
               break;

            case 3:
               spVector->setStateType(Equations::SphereExactStateIds::BENCHVELC2);
               break;
         }

         // Add magnetic initial state generator
         spVector = spGen->addVectorEquation<Equations::SphereExactVectorState>();
         spVector->setIdentity(PhysicalNames::MAGNETIC);
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
               tSH.clear(); 
               // Poloidal
               ptSH = tSH.insert(std::make_pair(std::make_pair(2,0), std::map<int,MHDComplex>()));
               ptSH.first->second.insert(std::make_pair(7, MHDComplex(1.0)));
               spVector->setHarmonicOptions(FieldComponents::Spectral::POL, tSH);
               break;

            case 3:
               spVector->setStateType(Equations::SphereExactStateIds::BENCHMAGC2);
               break;
         }

      // Generate random spectrum
      } else
      {
         // Shared pointer to random initial state equation
         Equations::SharedRandomVectorState spVector;

         // Add velocity random initial state generator 
         spVector = spGen->addVectorEquation<Equations::RandomVectorState>();
         spVector->setIdentity(PhysicalNames::VELOCITY);
         spVector->setSpectrum(FieldComponents::Spectral::TOR, -1e-4, 1e-4, 1e4, 1e4, 1e4);
         spVector->setSpectrum(FieldComponents::Spectral::POL, -1e-4, 1e-4, 1e4, 1e4, 1e4);

         // Add magnetic random initial state generator 
         spVector = spGen->addVectorEquation<Equations::RandomVectorState>();
         spVector->setIdentity(PhysicalNames::MAGNETIC);
         spVector->setSpectrum(FieldComponents::Spectral::TOR, -1e-4, 1e-4, 1e4, 1e4, 1e4);
         spVector->setSpectrum(FieldComponents::Spectral::POL, -1e-4, 1e-4, 1e4, 1e4, 1e4);
      }

      // Add output file
      IoVariable::SharedStateFileWriter spOut(new IoVariable::StateFileWriter(SchemeType::type(), SchemeType::isRegular()));
      spOut->expect(PhysicalNames::VELOCITY);
      spOut->expect(PhysicalNames::MAGNETIC);
      spGen->addHdf5OutputFile(spOut);
   }

   void BoussinesqPrecessionDynamoSphereStdModel::addVisualizers(SharedVisualizationGenerator spVis)
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

      // Add output file
      IoVariable::SharedVisualizationFileWriter spOut(new IoVariable::VisualizationFileWriter(SchemeType::type()));
      spOut->expect(PhysicalNames::VELOCITY);
      spOut->expect(PhysicalNames::VELOCITYZ);
      spOut->expect(PhysicalNames::VORTICITYZ);
      spOut->expect(PhysicalNames::MAGNETIC);
      spVis->addHdf5OutputFile(spOut);
   }

   void BoussinesqPrecessionDynamoSphereStdModel::setVisualizationState(SharedVisualizationGenerator spVis)
   {
      // Create and add initial state file to IO
      IoVariable::SharedStateFileReader spIn(new IoVariable::StateFileReader("4Visu", SchemeType::type(), SchemeType::isRegular()));

      // Set expected fields
      spIn->expect(PhysicalNames::VELOCITY);
      spIn->expect(PhysicalNames::MAGNETIC);

      // Set simulation state
      spVis->setInitialState(spIn);
   }

   void BoussinesqPrecessionDynamoSphereStdModel::addAsciiOutputFiles(SharedSimulation spSim)
   {
      // Create kinetic energy writer
      IoVariable::SharedSphereTorPolEnergyWriter spVector(new IoVariable::SphereTorPolEnergyWriter("kinetic", SchemeType::type()));
      spVector->expect(PhysicalNames::VELOCITY);
      spSim->addAsciiOutputFile(spVector);

      // Create magnetic energy writer
      spVector = IoVariable::SharedSphereTorPolEnergyWriter(new IoVariable::SphereTorPolEnergyWriter("magnetic", SchemeType::type()));
      spVector->expect(PhysicalNames::MAGNETIC);
      spSim->addAsciiOutputFile(spVector);
   }

   void BoussinesqPrecessionDynamoSphereStdModel::addHdf5OutputFiles(SharedSimulation spSim)
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

   void BoussinesqPrecessionDynamoSphereStdModel::addStatsOutputFiles(SharedSimulation spSim)
   {
   }

   void BoussinesqPrecessionDynamoSphereStdModel::setInitialState(SharedSimulation spSim)
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
