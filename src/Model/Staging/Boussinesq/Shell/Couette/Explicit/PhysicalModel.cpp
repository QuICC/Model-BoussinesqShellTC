/** 
 * @file PhysicalModel.cpp
 * @brief Source of the Boussinesq spherical Couette in a spherical shell (Toroidal/Poloidal formulation) without coupled solve (standard implementation)
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
#include MAKE_STR( QUICC_MODEL_PATH/Boussinesq/Shell/Couette/Explicit/PhysicalModel.hpp )

// Project includes
//
#include MAKE_STR( QUICC_MODEL_PATH/Boussinesq/Shell/Couette/Momentum.hpp )
#include "Enums/FieldIds.hpp"
#include "IoVariable/StateFileReader.hpp"
#include "IoVariable/StateFileWriter.hpp"
#include "IoVariable/VisualizationFileWriter.hpp"
#include "IoTools/IdToHuman.hpp"
#include "IoVariable/ShellTorPolEnergyWriter.hpp"
#include "IoVariable/ShellTorPolProbeWriter.hpp"
#include "IoVariable/ShellTorPolEnergySpectraWriter.hpp"
#include "IoVariable/ShellTorPolTorqueWriter.hpp"
#include "IoVariable/ShellTorPolUniformVorticityWriter.hpp"
#include "IoVariable/ShellTorPolDissipationWriter.hpp"
#include "IoVariable/ShellTorPolDissipationSpectraWriter.hpp"
#include "Generator/States/RandomVectorState.hpp"
#include "Generator/States/ShellExactStateIds.hpp"
#include "Generator/States/ShellExactVectorState.hpp"
#include "Generator/Visualizers/VectorFieldVisualizer.hpp"
#include "Generator/Visualizers/VectorFieldTrivialVisualizer.hpp"
#include "Generator/Visualizers/SphericalVerticalFieldVisualizer.hpp"
#include "Generator/Visualizers/SphericalRadialCylindricalFieldVisualizer.hpp"
#include "Model/PhysicalModelBase.hpp"

namespace QuICC {

namespace Model {

namespace Boussinesq {

namespace Shell {

namespace Couette {

namespace Explicit {

   const std::string PhysicalModel::PYMODULE = "boussinesq_couetteshell_std";

   const std::string PhysicalModel::PYCLASS = "BoussinesqCouetteShellStd";

   void PhysicalModel::addEquations(SharedSimulation spSim)
   {
      // Add Navier-Stokes equation
      spSim->addVectorEquation<Equations::Boussinesq::Shell::Couette::Momentum>();
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
            	spVector->setStateType(Equations::ShellExactStateIds::BENCHOMEGAX);
            	break;

            case 5:
            	spVector->setStateType(Equations::ShellExactStateIds::BENCHOMEGAY);
            	break;

            case 6:
            	spVector->setStateType(Equations::ShellExactStateIds::BENCHOMEGAZ);
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
      Equations::SharedVectorFieldTrivialVisualizer spVTrivial;
      Equations::SharedSphericalVerticalFieldVisualizer spVertical;
      Equations::SharedSphericalRadialCylindricalFieldVisualizer spCylindrical;

      // Add velocity field visualization
      spVector = spVis->addVectorEquation<Equations::VectorFieldVisualizer>();
      spVector->setFields(true, false, false);
      spVector->setIdentity(PhysicalNames::VELOCITY);

      // Add zonal velocity field visualization
      spVTrivial = spVis->addVectorEquation<Equations::VectorFieldTrivialVisualizer>();
      spVTrivial->setFields(true, false, false);
      spVTrivial->setIdentity(PhysicalNames::ZONAL_VELOCITY);

      // Add nonzonal velocity field visualization
      spVTrivial = spVis->addVectorEquation<Equations::VectorFieldTrivialVisualizer>();
      spVTrivial->setFields(true, false, false);
      spVTrivial->setIdentity(PhysicalNames::NONZONAL_VELOCITY);

      // Add vertical velocity visualization
      spVertical = spVis->addScalarEquation<Equations::SphericalVerticalFieldVisualizer>();
      spVertical->setFieldType(FieldType::VECTOR);
      spVertical->setIdentity(PhysicalNames::VELOCITYZ, PhysicalNames::VELOCITY);


      // Add cylindrical radius velocity visualization
      spCylindrical = spVis->addScalarEquation<Equations::SphericalRadialCylindricalFieldVisualizer>();
      spCylindrical->setFieldType(FieldType::VECTOR);
      spCylindrical->setIdentity(PhysicalNames::VELOCITYS, PhysicalNames::VELOCITY);

      // Add vertical vorticity visualization
      spVertical = spVis->addScalarEquation<Equations::SphericalVerticalFieldVisualizer>();
      spVertical->setFieldType(FieldType::CURL);
      spVertical->setIdentity(PhysicalNames::VORTICITYZ, PhysicalNames::VELOCITY);

      // Add horizontal vorticity visualization
      spCylindrical = spVis->addScalarEquation<Equations::SphericalRadialCylindricalFieldVisualizer>();
      spCylindrical->setFieldType(FieldType::CURL);
      spCylindrical->setIdentity(PhysicalNames::VORTICITYS, PhysicalNames::VELOCITY);


      // Add output file
      IoVariable::SharedVisualizationFileWriter spOut(new IoVariable::VisualizationFileWriter(SchemeType::type()));
      spOut->expect(PhysicalNames::VELOCITY);
      spOut->expect(PhysicalNames::ZONAL_VELOCITY);
      spOut->expect(PhysicalNames::NONZONAL_VELOCITY);
      spOut->expect(PhysicalNames::VELOCITYZ);
      spOut->expect(PhysicalNames::VELOCITYS);
      spOut->expect(PhysicalNames::VORTICITYZ);
      spOut->expect(PhysicalNames::VORTICITYS);
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
      IoVariable::SharedShellTorPolEnergyWriter spVector(new IoVariable::ShellTorPolEnergyWriter("kinetic", SchemeType::type()));
      spVector->expect(PhysicalNames::VELOCITY);
      spSim->addAsciiOutputFile(spVector);


      // Create probes
      Matrix mProbes(6,3);

      // the r coordinate is in x\in[-1,1] space
      // the theta coordinate is the cos of colatitude
      // the phi coordinate is azimuthal in radiants
      mProbes << 0.0, 0.5, 0.0,
    		  0.0, -0.5, Math::PI/3.,
			  0.0, 0.5, Math::PI*2./3.,
			  0.0, -0.5, Math::PI,
			  0.0, 0.5, Math::PI*4./3.,
			  0.0, -0.5, Math::PI*5./3.;

      // Create probe field writer
      IoVariable::SharedShellTorPolProbeWriter spVector2(new  IoVariable::ShellTorPolProbeWriter("velocity_probe", SchemeType::type(), mProbes));
      spVector2->expect(PhysicalNames::VELOCITY);
      spSim->addAsciiOutputFile(spVector2);


      // Create kinetic energy spectral writer
      IoVariable::SharedShellTorPolEnergySpectraWriter spVector3(new IoVariable::ShellTorPolEnergySpectraWriter("spectrum_kinetic", SchemeType::type()));
      spVector3->expect(PhysicalNames::VELOCITY);
      spSim->addAsciiOutputFile(spVector3);

      // Create torque writer
      IoVariable::SharedShellTorPolTorqueWriter spVector4(new IoVariable::ShellTorPolTorqueWriter("torque", SchemeType::type()));
      spVector4->expect(PhysicalNames::VELOCITY);
      spSim->addAsciiOutputFile(spVector4);

      // Create average vorticity writer
      IoVariable::SharedShellTorPolUniformVorticityWriter spVector5(new IoVariable::ShellTorPolUniformVorticityWriter("vorticity", SchemeType::type()));
      spVector5->expect(PhysicalNames::VELOCITY);
      spSim->addAsciiOutputFile(spVector5);

      // Create  kinetic dissipation writer
      IoVariable::SharedShellTorPolDissipationWriter spVector6(new IoVariable::ShellTorPolDissipationWriter("kinetic", SchemeType::type()));
      spVector6->expect(PhysicalNames::VELOCITY);
      spSim->addAsciiOutputFile(spVector6);

      // Create kinetic dissipation spectral writer
      IoVariable::SharedShellTorPolDissipationSpectraWriter spVector7(new IoVariable::ShellTorPolDissipationSpectraWriter("spectrum_kinetic", SchemeType::type()));
      spVector7->expect(PhysicalNames::VELOCITY);
      spSim->addAsciiOutputFile(spVector7);

   }

   void PhysicalModel::addHdf5OutputFiles(SharedSimulation spSim)
   {
      // Field IDs iterator
      std::vector<PhysicalNames::Id>::const_iterator it;
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
