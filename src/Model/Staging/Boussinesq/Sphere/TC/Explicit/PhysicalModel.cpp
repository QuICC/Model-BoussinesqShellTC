/** 
 * @file PhysicalModel.cpp
 * @brief Source of the Boussinesq thermal convection in a sphere (Toroidal/Poloidal formulation) without coupled solve (standard implementation)
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
#include MAKE_STR( QUICC_MODEL_PATH/Boussinesq/Sphere/TC/Explicit/PhysicalModel.hpp )

// Project includes
//
#include MAKE_STR( QUICC_MODEL_PATH/Boussinesq/Sphere/TC/Transport.hpp )
#include MAKE_STR( QUICC_MODEL_PATH/Boussinesq/Sphere/TC/Momentum.hpp )
#include "Enums/FieldIds.hpp"
#include "IoVariable/StateFileReader.hpp"
#include "IoVariable/StateFileWriter.hpp"
#include "IoVariable/VisualizationFileWriter.hpp"
#include "IoTools/IdToHuman.hpp"
#include "IoVariable/SphereAngularMomentumWriter.hpp"
#include "IoVariable/SphereNusseltWriter.hpp"
#include "IoVariable/SphereScalarEnergyWriter.hpp"
#include "IoVariable/SphereScalarLSpectrumWriter.hpp"
#include "IoVariable/SphereScalarMSpectrumWriter.hpp"
#include "IoVariable/SphereTorPolEnergyWriter.hpp"
#include "IoVariable/SphereTorPolLSpectrumWriter.hpp"
#include "IoVariable/SphereTorPolMSpectrumWriter.hpp"
#include "IoVariable/SphereTorPolEnstrophyWriter.hpp"
#include "IoVariable/SphereTorPolEnstrophyLSpectrumWriter.hpp"
#include "IoVariable/SphereTorPolEnstrophyMSpectrumWriter.hpp"
#include "Generator/States/RandomScalarState.hpp"
#include "Generator/States/RandomVectorState.hpp"
#include "Generator/States/SphereExactStateIds.hpp"
#include "Generator/States/SphereExactScalarState.hpp"
#include "Generator/States/SphereExactVectorState.hpp"
#include "Generator/Visualizers/ScalarFieldVisualizer.hpp"
#include "Generator/Visualizers/VectorFieldVisualizer.hpp"
#include "Model/PhysicalModelBase.hpp"


namespace QuICC {

namespace Model {

namespace Boussinesq {

namespace Sphere {

namespace TC {

namespace Explicit {

   const std::string PhysicalModel::PYMODULE = "boussinesq_tcsphere_std";

   const std::string PhysicalModel::PYCLASS = "BoussinesqTCSphereStd";

   void PhysicalModel::addEquations(SharedSimulation spSim)
   {
      // Add transport equation
      spSim->addScalarEquation<Equations::Boussinesq::Sphere::TC::Transport>();
                                                              
      // Add Navier-Stokes equation                           
      spSim->addVectorEquation<Equations::Boussinesq::Sphere::TC::Momentum>();
   }

   void PhysicalModel::addStates(SharedStateGenerator spGen)
   {
      // Generate "exact" solutions (trigonometric or monomial)
      if(true)
      {
         // Shared pointer to equation
         Equations::SharedSphereExactScalarState spScalar;
         Equations::SharedSphereExactVectorState spVector;

         Equations::SphereExactVectorState::HarmonicModeType tSH;
         std::pair<Equations::SphereExactVectorState::HarmonicModeType::iterator,bool> ptSH; 

         // Add temperature initial state generator
         spScalar = spGen->addScalarEquation<Equations::SphereExactScalarState>();
         spScalar->setIdentity(PhysicalNames::TEMPERATURE);
         switch(0)
         {
            case 0:
               spScalar->setSpectralType(Equations::SphereExactStateIds::HARMONIC);
               tSH.clear(); 
               ptSH = tSH.insert(std::make_pair(std::make_pair(3,3), std::map<int,MHDComplex>()));
               ptSH.first->second.insert(std::make_pair(7, MHDComplex(1.0,2.0)));
               spScalar->setHarmonicOptions(tSH);
               break;
         }

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
	       //Exact
	       spVector->setStateType(Equations::SphereExactStateIds::VALIDATION_ENSTROPHY2);
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
         spScalar->setSpectrum(-1e-2, 1e-2, 1e4, 1e4, 1e4);
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
      IoVariable::SharedSphereScalarEnergyWriter spTemp(new IoVariable::SphereScalarEnergyWriter("temperature", SchemeType::type()));
      spTemp->expect(PhysicalNames::TEMPERATURE);
      spSim->addAsciiOutputFile(spTemp);

      // Create temperature L energy spectrum writer
      IoVariable::SharedSphereScalarLSpectrumWriter spTempL(new IoVariable::SphereScalarLSpectrumWriter("temperature", SchemeType::type()));
      spTempL->expect(PhysicalNames::TEMPERATURE);
      //spTempL->numberOutput();
      spSim->addAsciiOutputFile(spTempL);

      // Create temperature M energy spectrum writer
      IoVariable::SharedSphereScalarMSpectrumWriter spTempM(new IoVariable::SphereScalarMSpectrumWriter("temperature", SchemeType::type()));
      spTempM->expect(PhysicalNames::TEMPERATURE);
      //spTempM->numberOutput();
      spSim->addAsciiOutputFile(spTempM);


      // Create Nusselt number writer
      IoVariable::SharedSphereNusseltWriter spNusselt(new IoVariable::SphereNusseltWriter("", SchemeType::type()));
      spNusselt->expect(PhysicalNames::TEMPERATURE);
      spSim->addAsciiOutputFile(spNusselt);


      // Create kinetic energy writer
      IoVariable::SharedSphereTorPolEnergyWriter spKinetic(new IoVariable::SphereTorPolEnergyWriter("kinetic", SchemeType::type()));
      spKinetic->expect(PhysicalNames::VELOCITY);
      spSim->addAsciiOutputFile(spKinetic);

      // Create kinetic L energy spectrum writer
      IoVariable::SharedSphereTorPolLSpectrumWriter spKineticL(new IoVariable::SphereTorPolLSpectrumWriter("kinetic", SchemeType::type()));
      spKineticL->expect(PhysicalNames::VELOCITY);
      //spKineticL->numberOutput();
      spSim->addAsciiOutputFile(spKineticL);

      // Create kinetic M energy spectrum writer
      IoVariable::SharedSphereTorPolMSpectrumWriter spKineticM(new IoVariable::SphereTorPolMSpectrumWriter("kinetic", SchemeType::type()));
      spKineticM->expect(PhysicalNames::VELOCITY);
      //spKineticM->numberOutput();
      spSim->addAsciiOutputFile(spKineticM);
  
      // Create enstrophy writer
      IoVariable::SharedSphereTorPolEnstrophyWriter spEnstrophy(new IoVariable::SphereTorPolEnstrophyWriter("kinetic", SchemeType::type()));
      spEnstrophy->expect(PhysicalNames::VELOCITY);
      spSim->addAsciiOutputFile(spEnstrophy);

      // Create kinetic enstrophy L spectrum writer
      IoVariable::SharedSphereTorPolEnstrophyLSpectrumWriter spEnstrophyL(new IoVariable::SphereTorPolEnstrophyLSpectrumWriter("kinetic_enstrophy", SchemeType::type()));
      spEnstrophyL->expect(PhysicalNames::VELOCITY);
      //spEnstrophyL->numberOutput();
      spSim->addAsciiOutputFile(spEnstrophyL);

      // Create kinetic enstrophy M spectrum writer
      IoVariable::SharedSphereTorPolEnstrophyMSpectrumWriter spEnstrophyM(new IoVariable::SphereTorPolEnstrophyMSpectrumWriter("kinetic_enstrophy", SchemeType::type()));
      spEnstrophyM->expect(PhysicalNames::VELOCITY);
      //spEnstrophyM->numberOutput();
      spSim->addAsciiOutputFile(spEnstrophyM);

#if 1
      // Create angular momentum writer
      IoVariable::SharedSphereAngularMomentumWriter spAngMom(new IoVariable::SphereAngularMomentumWriter("", SchemeType::type()));
      spAngMom->expect(PhysicalNames::VELOCITY);
      spSim->addAsciiOutputFile(spAngMom);
#endif
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
