/** 
 * @file BoussinesqDynamoShellStdModel.cpp
 * @brief Source of the Boussinesq thermal convection dynamo in a spherical shell (Toroidal/Poloidal formulation) without coupled solve (standard implementation)
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
#include "PhysicalModels/BoussinesqDynamoShellStdModel.hpp"

// Project includes
//
#include "Enums/FieldIds.hpp"
#include "IoVariable/StateFileReader.hpp"
#include "IoVariable/StateFileWriter.hpp"
#include "IoVariable/VisualizationFileWriter.hpp"
#include "IoTools/IdToHuman.hpp"
#include "Equations/Shell/Boussinesq/BoussinesqDynamoShellTransport.hpp"
#include "Equations/Shell/Boussinesq/BoussinesqDynamoShellMomentum.hpp"
#include "Equations/Shell/Boussinesq/BoussinesqDynamoShellInduction.hpp"
#include "Generator/States/RandomScalarState.hpp"
#include "Generator/States/RandomVectorState.hpp"
#include "Generator/States/ShellExactStateIds.hpp"
#include "Generator/States/ShellExactScalarState.hpp"
#include "Generator/States/ShellExactVectorState.hpp"
#include "Generator/Visualizers/ScalarFieldVisualizer.hpp"
#include "Generator/Visualizers/VectorFieldVisualizer.hpp"
#include "PhysicalModels/PhysicalModelBase.hpp"

namespace GeoMHDiSCC {

   const std::string BoussinesqDynamoShellStdModel::PYMODULE = "boussinesq_dynamoshell_std";

   const std::string BoussinesqDynamoShellStdModel::PYCLASS = "BoussinesqDynamoShellStd";

   void BoussinesqDynamoShellStdModel::addEquations(SharedSimulation spSim)
   {
      // Add transport equation
      spSim->addScalarEquation<Equations::BoussinesqDynamoShellTransport>();
      
      // Add Navier-Stokes equation
      spSim->addVectorEquation<Equations::BoussinesqDynamoShellMomentum>();
      
      // Add induction equation
      spSim->addVectorEquation<Equations::BoussinesqDynamoShellInduction>();
   }

   void BoussinesqDynamoShellStdModel::addStates(SharedStateGenerator spGen)
   {
      // Generate "exact" solutions (trigonometric or monomial)
      if(true)
      {
         // Shared pointer to equation
         Equations::SharedShellExactScalarState spScalar;
         Equations::SharedShellExactVectorState spVector;

         // Add temperature initial state generator
         spScalar = spGen->addScalarEquation<Equations::ShellExactScalarState>();
         spScalar->setIdentity(PhysicalNames::TEMPERATURE);
         spScalar->setStateType(Equations::ShellExactStateIds::HARMONIC);
         std::vector<std::tr1::tuple<int,int,MHDComplex> > tSH;
         tSH.push_back(std::tr1::make_tuple(5,5,MHDComplex(1,1)));
         spScalar->setHarmonicOptions(tSH);

         // Add velocity initial state generator
         spVector = spGen->addVectorEquation<Equations::ShellExactVectorState>();
         spVector->setIdentity(PhysicalNames::VELOCITY);
         switch(1)
         {
            case 0:
               spVector->setStateType(FieldComponents::Physical::R, Equations::ShellExactStateIds::TORPOLT11P11);
               spVector->setStateType(FieldComponents::Physical::THETA, Equations::ShellExactStateIds::TORPOLT11P11);
               spVector->setStateType(FieldComponents::Physical::PHI, Equations::ShellExactStateIds::TORPOLT11P11);
               break;

            case 1:
               spVector->setStateType(FieldComponents::Physical::R, Equations::ShellExactStateIds::TORPOLT54P43);
               spVector->setStateType(FieldComponents::Physical::THETA, Equations::ShellExactStateIds::TORPOLT54P43);
               spVector->setStateType(FieldComponents::Physical::PHI, Equations::ShellExactStateIds::TORPOLT54P43);
               break;

            case 2:
               spVector->setStateType(FieldComponents::Physical::R, Equations::ShellExactStateIds::HARMONIC);
               tSH.clear(); 
               tSH.push_back(std::tr1::make_tuple(3,1,MHDComplex(1,0)));
               spVector->setHarmonicOptions(FieldComponents::Physical::R, tSH);
               spVector->setStateType(FieldComponents::Physical::THETA, Equations::ShellExactStateIds::HARMONIC);
               tSH.clear(); 
               tSH.push_back(std::tr1::make_tuple(1,1,MHDComplex(1,0)));
               spVector->setHarmonicOptions(FieldComponents::Physical::THETA, tSH);
               spVector->setStateType(FieldComponents::Physical::PHI, Equations::ShellExactStateIds::HARMONIC);
               tSH.clear(); 
               tSH.push_back(std::tr1::make_tuple(2,1,MHDComplex(1,0)));
               spVector->setHarmonicOptions(FieldComponents::Physical::PHI, tSH);
         }

         // Add magnetic initial state generator
         spVector = spGen->addVectorEquation<Equations::ShellExactVectorState>();
         spVector->setIdentity(PhysicalNames::MAGNETIC);
         switch(1)
         {
            case 0:
               spVector->setStateType(FieldComponents::Physical::R, Equations::ShellExactStateIds::TORPOLT11P11);
               spVector->setStateType(FieldComponents::Physical::THETA, Equations::ShellExactStateIds::TORPOLT11P11);
               spVector->setStateType(FieldComponents::Physical::PHI, Equations::ShellExactStateIds::TORPOLT11P11);
               break;

            case 1:
               spVector->setStateType(FieldComponents::Physical::R, Equations::ShellExactStateIds::TORPOLT54P43);
               spVector->setStateType(FieldComponents::Physical::THETA, Equations::ShellExactStateIds::TORPOLT54P43);
               spVector->setStateType(FieldComponents::Physical::PHI, Equations::ShellExactStateIds::TORPOLT54P43);
               break;

            case 2:
               spVector->setStateType(FieldComponents::Physical::R, Equations::ShellExactStateIds::HARMONIC);
               tSH.clear(); 
               tSH.push_back(std::tr1::make_tuple(3,1,MHDComplex(1,0)));
               spVector->setHarmonicOptions(FieldComponents::Physical::R, tSH);
               spVector->setStateType(FieldComponents::Physical::THETA, Equations::ShellExactStateIds::HARMONIC);
               tSH.clear(); 
               tSH.push_back(std::tr1::make_tuple(1,1,MHDComplex(1,0)));
               spVector->setHarmonicOptions(FieldComponents::Physical::THETA, tSH);
               spVector->setStateType(FieldComponents::Physical::PHI, Equations::ShellExactStateIds::HARMONIC);
               tSH.clear(); 
               tSH.push_back(std::tr1::make_tuple(2,1,MHDComplex(1,0)));
               spVector->setHarmonicOptions(FieldComponents::Physical::PHI, tSH);
         }

      // Generate random spectrum
      } else
      {
         // Shared pointer to random initial state equation
         Equations::SharedRandomScalarState spScalar;
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

         // Add scalar random initial state generator
         spScalar = spGen->addScalarEquation<Equations::RandomScalarState>();
         spScalar->setIdentity(PhysicalNames::TEMPERATURE);
         spScalar->setSpectrum(-1e-3, 1e-3, 1e4, 1e4, 1e4);
      }

      // Add output file
      IoVariable::SharedStateFileWriter spOut(new IoVariable::StateFileWriter(SchemeType::type(), SchemeType::isRegular()));
      spOut->expect(PhysicalNames::TEMPERATURE);
      spOut->expect(PhysicalNames::VELOCITY);
      spOut->expect(PhysicalNames::MAGNETIC);
      spGen->addHdf5OutputFile(spOut);
   }

   void BoussinesqDynamoShellStdModel::addVisualizers(SharedVisualizationGenerator spVis)
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

      // Add magnetic field visualization
      spVector = spVis->addVectorEquation<Equations::VectorFieldVisualizer>();
      spVector->setFields(true, false, true);
      spVector->setIdentity(PhysicalNames::MAGNETIC);

      // Add output file
      IoVariable::SharedVisualizationFileWriter spOut(new IoVariable::VisualizationFileWriter(SchemeType::type()));
      spOut->expect(PhysicalNames::TEMPERATURE);
      spOut->expect(PhysicalNames::VELOCITY);
      spOut->expect(PhysicalNames::MAGNETIC);
      spVis->addHdf5OutputFile(spOut);
   }

   void BoussinesqDynamoShellStdModel::setVisualizationState(SharedVisualizationGenerator spVis)
   {
      // Create and add initial state file to IO
      IoVariable::SharedStateFileReader spIn(new IoVariable::StateFileReader("4Visu", SchemeType::type(), SchemeType::isRegular()));

      // Set expected fields
      spIn->expect(PhysicalNames::TEMPERATURE);
      spIn->expect(PhysicalNames::VELOCITY);
      spIn->expect(PhysicalNames::MAGNETIC);

      // Set simulation state
      spVis->setInitialState(spIn);
   }

   void BoussinesqDynamoShellStdModel::addAsciiOutputFiles(SharedSimulation spSim)
   {
      // Add ASCII output file
      //pSim->addOutputFile(AN_ASCIIFILE);
      
      // Add ASCII output file
      //pSim->addOutputFile(AN_ASCIIFILE);
   }

   void BoussinesqDynamoShellStdModel::addHdf5OutputFiles(SharedSimulation spSim)
   {
      // Field IDs iterator
      std::vector<GeoMHDiSCC::PhysicalNames::Id>::const_iterator  it;
      std::vector<GeoMHDiSCC::PhysicalNames::Id> ids = PhysicalModelBase::fieldIds();

      // Create and add state file to IO
      IoVariable::SharedStateFileWriter spState(new IoVariable::StateFileWriter(SchemeType::type(), SchemeType::isRegular()));
      for(it = ids.begin(); it != ids.end(); ++it)
      {
         spState->expect(*it);
      }
      spSim->addHdf5OutputFile(spState);
   }

   void BoussinesqDynamoShellStdModel::setInitialState(SharedSimulation spSim)
   {
      // Field IDs iterator
      std::vector<GeoMHDiSCC::PhysicalNames::Id>::const_iterator  it;
      std::vector<GeoMHDiSCC::PhysicalNames::Id> ids = PhysicalModelBase::fieldIds();

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
