/** 
 * @file PhysicalModel.cpp
 * @brief Source of the Boussinesq F-plane 3DQG physical model
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
#include MAKE_STR( QUICC_MODEL_PATH/Boussinesq/Plane/F3DQG/PhysicalModel.hpp )

// Project includes
//
#include MAKE_STR( QUICC_MODEL_PATH/Boussinesq/Plane/F3DQG/Streamfunction.hpp )
#include MAKE_STR( QUICC_MODEL_PATH/Boussinesq/Plane/F3DQG/VelocityX.hpp )
#include MAKE_STR( QUICC_MODEL_PATH/Boussinesq/Plane/F3DQG/VelocityY.hpp )
#include MAKE_STR( QUICC_MODEL_PATH/Boussinesq/Plane/F3DQG/VelocityZ.hpp )
#include MAKE_STR( QUICC_MODEL_PATH/Boussinesq/Plane/F3DQG/Transport.hpp )
#include MAKE_STR( QUICC_MODEL_PATH/Boussinesq/Plane/F3DQG/VorticityZ.hpp )
#include MAKE_STR( QUICC_MODEL_PATH/Boussinesq/Plane/F3DQG/MeanHeat.hpp )
#include MAKE_STR( QUICC_MODEL_PATH/Boussinesq/Plane/F3DQG/DissTh.hpp )
#include MAKE_STR( QUICC_MODEL_PATH/Boussinesq/Plane/F3DQG/DissV.hpp )
#include "Enums/FieldIds.hpp"
#include "IoVariable/StateFileReader.hpp"
#include "IoVariable/StateFileWriter.hpp"
#include "IoVariable/VisualizationFileWriter.hpp"
#include "IoTools/IdToHuman.hpp"
#include "IoVariable/Cartesian1DNusseltZWriter.hpp"
#include "IoVariable/Cartesian1DScalarEnergyWriter.hpp"
#include "IoVariable/Cartesian1DStreamEnergyWriter.hpp"
#include "IoVariable/Cartesian1DMagneticEnergyWriter.hpp"
#include "IoVariable/Cartesian1DFluctuatingMagneticEnergyWriter.hpp"
#include "IoVariable/Cartesian1DKineticCartesianWriter.hpp"
#include "IoStats/Cartesian1DScalarAvgWriter.hpp"
#include "IoStats/Cartesian1DScalarRMSWriter.hpp"
#include "IoStats/Cartesian1DScalarSkewWriter.hpp"
#include "IoStats/Cartesian1DScalarKurtWriter.hpp"
#include "Generator/States/RandomScalarState.hpp"
#include "Generator/States/CartesianExactScalarState.hpp"
#include "Generator/Visualizers/ScalarFieldVisualizer.hpp"
#include "Model/PhysicalModelBase.hpp"

namespace QuICC {

namespace Model {

namespace Boussinesq {

namespace Plane {

namespace F3DQG {

   const std::string PhysicalModel::PYMODULE = "boussinesq_fplane3dqg";

   const std::string PhysicalModel::PYCLASS = "BoussinesqFPlane3DQG";

   void PhysicalModel::addEquations(SharedSimulation spSim)
   {
      // Add upright streamfunction equation
      spSim->addScalarEquation<Equations::Boussinesq::Plane::F3DQG::Streamfunction>();
      
      // Add upright vertical velocity equation
      spSim->addScalarEquation<Equations::Boussinesq::Plane::F3DQG::VelocityZ>();
      
      // Add upright transport equation
      spSim->addScalarEquation<Equations::Boussinesq::Plane::F3DQG::Transport>();

      // Add velocity x
      spSim->addScalarEquation<Equations::Boussinesq::Plane::F3DQG::VelocityX>();

      // Add velocity y
      spSim->addScalarEquation<Equations::Boussinesq::Plane::F3DQG::VelocityY>();

      // Add vertical vorticity equation
      spSim->addScalarEquation<Equations::Boussinesq::Plane::F3DQG::VorticityZ>();

      // Add mean heat computation
      spSim->addScalarEquation<Equations::Boussinesq::Plane::F3DQG::MeanHeat>();

      // Add DissV computation
      spSim->addScalarEquation<Equations::Boussinesq::Plane::F3DQG::DissV>();

      // Add DissTh computation
      spSim->addScalarEquation<Equations::Boussinesq::Plane::F3DQG::DissTh>();
   }

   void PhysicalModel::addStates(SharedStateGenerator spGen)
   {
      // Generate "exact" solutions (trigonometric or monomial)
      if(false)
      {
         // Shared pointer to equation
         Equations::SharedCartesianExactScalarState spExact;

         // Add transport initial state generation equation
         spExact = spGen->addScalarEquation<Equations::CartesianExactScalarState>();
         spExact->setIdentity(PhysicalNames::TEMPERATURE);
         spExact->setStateType(Equations::CartesianExactStateIds::POLYCOSCOS);
         spExact->setModeOptions(1e0, 0.0, 1e0, 0.0, 1e0, 0.0);

         // Add streamfunction initial state generation equation
         spExact = spGen->addScalarEquation<Equations::CartesianExactScalarState>();
         spExact->setIdentity(PhysicalNames::STREAMFUNCTION);
         spExact->setStateType(Equations::CartesianExactStateIds::POLYCOSCOS);
         spExact->setModeOptions(2e0, 0.0, 1e0, 0.0, 1e0, 0.0);

         // Add vertical velocity initial state generation equation
         spExact = spGen->addScalarEquation<Equations::CartesianExactScalarState>();
         spExact->setIdentity(PhysicalNames::VELOCITYZ);
         spExact->setStateType(Equations::CartesianExactStateIds::POLYCOSCOS);
         spExact->setModeOptions(3e0, 0.0, 1e0, 0.0, 1e0, 0.0);

         // Add vertical velocity initial state generation equation
         spExact = spGen->addScalarEquation<Equations::CartesianExactScalarState>();
         spExact->setIdentity(PhysicalNames::DZ_MEANTEMPERATURE);
         spExact->setStateType(Equations::CartesianExactStateIds::POLYCOSCOS);
         spExact->setModeOptions(-1e0, 0.0, 1e0, 0.0, 1e0, 0.0);

      // Generate random spectrum
      } else
      {
         // Shared pointer to random initial state equation
         Equations::SharedRandomScalarState spRand;
         Equations::SharedCartesianExactScalarState spExact;

         // Add transport initial state generation equation
         spRand = spGen->addScalarEquation<Equations::RandomScalarState>();
         spRand->setIdentity(PhysicalNames::TEMPERATURE);
         spRand->setSpectrum(-1e-2, 1e-2, 1e4, 1e4, 1e4);

         // Add streamfunction initial state generation equation
         spRand = spGen->addScalarEquation<Equations::RandomScalarState>();
         spRand->setIdentity(PhysicalNames::STREAMFUNCTION);
         spRand->setSpectrum(-1e-2, 1e-2, 1e4, 1e4, 1e4);

         // Add vertical velocity initial state generation equation
         spRand = spGen->addScalarEquation<Equations::RandomScalarState>();
         spRand->setIdentity(PhysicalNames::VELOCITYZ);
         spRand->setSpectrum(-1e-2, 1e-2, 1e4, 1e4, 1e4);

         // Add vertical velocity initial state generation equation
         spExact = spGen->addScalarEquation<Equations::CartesianExactScalarState>();
         spExact->setIdentity(PhysicalNames::DZ_MEANTEMPERATURE);
         spExact->setStateType(Equations::CartesianExactStateIds::POLYCOSCOS);
         spExact->setModeOptions(0e0, 0.0, 1e0, 0.0, 1e0, 0.0);
      }

      // Add output file
      IoVariable::SharedStateFileWriter spOut(new IoVariable::StateFileWriter(SchemeType::type(), SchemeType::isRegular()));
      spOut->expect(PhysicalNames::TEMPERATURE);
      spOut->expect(PhysicalNames::STREAMFUNCTION);
      spOut->expect(PhysicalNames::VELOCITYZ);
      spOut->expect(PhysicalNames::DZ_MEANTEMPERATURE);
      spGen->addHdf5OutputFile(spOut);
   }

   void PhysicalModel::addVisualizers(SharedVisualizationGenerator spVis)
   {
      // Shared pointer to basic field visualizer
      Equations::SharedScalarFieldVisualizer spField;

      // Add transport field visualization
      spField = spVis->addScalarEquation<Equations::ScalarFieldVisualizer>();
      spField->setFields(true, false, false);
      spField->setIdentity(PhysicalNames::TEMPERATURE);
      
      // Add streamfunction field visualization
      spField = spVis->addScalarEquation<Equations::ScalarFieldVisualizer>();
      spField->setFields(true, true, false);
      spField->setIdentity(PhysicalNames::STREAMFUNCTION);

      // Add VelocityX to profile visualisation
      spField = spVis->addScalarEquation<Equations::ScalarFieldVisualizer>();
      spField->setFields(true, false);
      spField->setIdentity(PhysicalNames::VELOCITYX);

      // Add VelocityY to profile visualisation
      spField = spVis->addScalarEquation<Equations::ScalarFieldVisualizer>();
      spField->setFields(true, false);
      spField->setIdentity(PhysicalNames::VELOCITYY);
      
      // Add vertical velocity field visualization
      spField = spVis->addScalarEquation<Equations::ScalarFieldVisualizer>();
      spField->setFields(true, false, false);
      spField->setIdentity(PhysicalNames::VELOCITYZ);
      
      // Add background temperature profile visualization
      spField = spVis->addScalarEquation<Equations::ScalarFieldVisualizer>();
      spField->setFields(true, false, false);
      spField->setIdentity(PhysicalNames::DZ_MEANTEMPERATURE);

      // Add vorticity profile visualization
      spField = spVis->addScalarEquation<Equations::ScalarFieldVisualizer>();
      spField->setFields(true, false, false);
      spField->setIdentity(PhysicalNames::VORTICITYZ);

      // Add output file
      IoVariable::SharedVisualizationFileWriter spOut(new IoVariable::VisualizationFileWriter(SchemeType::type()));
      spOut->expect(PhysicalNames::TEMPERATURE);
      spOut->expect(PhysicalNames::STREAMFUNCTION);
      spOut->expect(PhysicalNames::VELOCITYX);
      spOut->expect(PhysicalNames::VELOCITYY);
      spOut->expect(PhysicalNames::VELOCITYZ);
      spOut->expect(PhysicalNames::DZ_MEANTEMPERATURE);
      spOut->expect(PhysicalNames::VORTICITYZ);
      spVis->addHdf5OutputFile(spOut);
   }

   void PhysicalModel::setVisualizationState(SharedVisualizationGenerator spVis)
   {
      // Create and add initial state file to IO
      IoVariable::SharedStateFileReader spIn(new IoVariable::StateFileReader("4Visu", SchemeType::type(), SchemeType::isRegular()));

      // Set expected fields
      spIn->expect(PhysicalNames::TEMPERATURE);
      spIn->expect(PhysicalNames::STREAMFUNCTION);
      spIn->expect(PhysicalNames::VELOCITYX);
      spIn->expect(PhysicalNames::VELOCITYY);
      spIn->expect(PhysicalNames::VELOCITYZ);
      spIn->expect(PhysicalNames::DZ_MEANTEMPERATURE);
      spIn->expect(PhysicalNames::VORTICITYZ);

      // Set simulation state
      spVis->setInitialState(spIn);
   }

   void PhysicalModel::addAsciiOutputFiles(SharedSimulation spSim)
   {
      // Create Nusselt number writer
      IoVariable::SharedCartesian1DNusseltZWriter spNusselt(new IoVariable::Cartesian1DNusseltZWriter(SchemeType::type()));
      spNusselt->expect(PhysicalNames::DZ_MEANTEMPERATURE);
      spSim->addAsciiOutputFile(spNusselt);

      // Create temperature energy writer
      IoVariable::SharedCartesian1DScalarEnergyWriter spTemp(new IoVariable::Cartesian1DScalarEnergyWriter("temperature", SchemeType::type()));
      spTemp->expect(PhysicalNames::TEMPERATURE);
      spSim->addAsciiOutputFile(spTemp);

      // Create kinetic energy writer
      IoVariable::SharedCartesian1DStreamEnergyWriter spStream(new IoVariable::Cartesian1DStreamEnergyWriter("kinetic", SchemeType::type()));
      spStream->expect(PhysicalNames::STREAMFUNCTION);
      spStream->expect(PhysicalNames::VELOCITYZ);
      spSim->addAsciiOutputFile(spStream);

      // Create cartesian kinetic energy writer
      // to be compared with kinetic_energy, created by calculating velocityx and velocityy
      // NB: It's mangetic energy per unit volume
      IoVariable::SharedCartesian1DKineticCartesianWriter spCartK(new IoVariable::Cartesian1DKineticCartesianWriter("cartesian_kinetic", SchemeType::type()));
      spCartK->expect(PhysicalNames::VELOCITYX);
      spCartK->expect(PhysicalNames::VELOCITYY);
      spCartK->expect(PhysicalNames::VELOCITYZ);
      spSim->addAsciiOutputFile(spCartK);
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

      // Add mean temperature to ouput file
      spState->expect(PhysicalNames::DZ_MEANTEMPERATURE);
      spState->expect(PhysicalNames::VORTICITYZ);
      spState->expect(PhysicalNames::VELOCITYX);
      spState->expect(PhysicalNames::VELOCITYY);

      spSim->addHdf5OutputFile(spState);
   }

   void PhysicalModel::addStatsOutputFiles(SharedSimulation spSim)
   {
      // Create several stats (Avg, RMS, Skew and Kurt) for various fields

      // temperature

      // Create Avg temperature writer
      IoStats::SharedCartesian1DScalarAvgWriter spAvg(new IoStats::Cartesian1DScalarAvgWriter("temperature",SchemeType::type()));
      spAvg->expect(PhysicalNames::TEMPERATURE);
      spSim->addStatsOutputFile(spAvg);

      // Create RMS temperature writer
       IoStats::SharedCartesian1DScalarRMSWriter spRMS(new IoStats::Cartesian1DScalarRMSWriter("temperature", spAvg,  SchemeType::type()));
      spRMS->expect(PhysicalNames::TEMPERATURE);
      spSim->addStatsOutputFile(spRMS);

     // Create skew temperature writer
       IoStats::SharedCartesian1DScalarSkewWriter spSkew(new IoStats::Cartesian1DScalarSkewWriter("temperature", spAvg, spRMS,  SchemeType::type()));
      spSkew->expect(PhysicalNames::TEMPERATURE);
      spSim->addStatsOutputFile(spSkew);

      // Create kurt temperature writer
      IoStats::SharedCartesian1DScalarKurtWriter spKurt(new IoStats::Cartesian1DScalarKurtWriter("temperature", spAvg, spRMS,  SchemeType::type()));
      spKurt->expect(PhysicalNames::TEMPERATURE);
      spSim->addStatsOutputFile(spKurt);

      // dz_meantemperature

      // Create Avg dz(mean temperature) writer
      IoStats::SharedCartesian1DScalarAvgWriter spAvgDZMT(new IoStats::Cartesian1DScalarAvgWriter("dz_meantemperature",SchemeType::type()));
      spAvgDZMT->expect(PhysicalNames::DZ_MEANTEMPERATURE);
      spSim->addStatsOutputFile(spAvgDZMT);

      // Create RMS dz(mean temperature) writer
      IoStats::SharedCartesian1DScalarRMSWriter spRMSDZMT(new IoStats::Cartesian1DScalarRMSWriter("dz_meantemperature", spAvgDZMT,  SchemeType::type()));
      spRMSDZMT->expect(PhysicalNames::DZ_MEANTEMPERATURE);
      spSim->addStatsOutputFile(spRMSDZMT);

     // Create skew dz(mean temperature) writer
       IoStats::SharedCartesian1DScalarSkewWriter spSkewDZMT(new IoStats::Cartesian1DScalarSkewWriter("dz_meantemperature", spAvgDZMT, spRMSDZMT,  SchemeType::type()));
      spSkewDZMT->expect(PhysicalNames::DZ_MEANTEMPERATURE);
      spSim->addStatsOutputFile(spSkewDZMT);

      // Create kurt dz(mean temperature) writer
      IoStats::SharedCartesian1DScalarKurtWriter spKurtDZMT(new IoStats::Cartesian1DScalarKurtWriter("dz_meantemperature", spAvgDZMT, spRMSDZMT,  SchemeType::type()));
      spKurtDZMT->expect(PhysicalNames::DZ_MEANTEMPERATURE);
      spSim->addStatsOutputFile(spKurtDZMT);

      // velocityz

      // Create Avg vertical velocity writer
      IoStats::SharedCartesian1DScalarAvgWriter spAvgW(new IoStats::Cartesian1DScalarAvgWriter("velocityz",SchemeType::type()));
      spAvgW->expect(PhysicalNames::VELOCITYZ);
      spSim->addStatsOutputFile(spAvgW);

      // Create RMS vertical velocity writer
      IoStats::SharedCartesian1DScalarRMSWriter spRMSW(new IoStats::Cartesian1DScalarRMSWriter("velocityz", spAvgW,  SchemeType::type()));
      spRMSW->expect(PhysicalNames::VELOCITYZ);
      spSim->addStatsOutputFile(spRMSW);

     // Create skew vertical velocity writer
       IoStats::SharedCartesian1DScalarSkewWriter spSkewW(new IoStats::Cartesian1DScalarSkewWriter("velocityz", spAvgW, spRMSW,  SchemeType::type()));
      spSkewW->expect(PhysicalNames::VELOCITYZ);
      spSim->addStatsOutputFile(spSkewW);

      // Create kurt vertical velocity writer
      IoStats::SharedCartesian1DScalarKurtWriter spKurtW(new IoStats::Cartesian1DScalarKurtWriter("velocityz", spAvgW, spRMSW,  SchemeType::type()));
      spKurtW->expect(PhysicalNames::VELOCITYZ);
      spSim->addStatsOutputFile(spKurtW);

      // vorticityz

      // Create Avg vertical vorticity writer
      IoStats::SharedCartesian1DScalarAvgWriter spAvgZ(new IoStats::Cartesian1DScalarAvgWriter("vorticityz",SchemeType::type()));
      spAvgZ->expect(PhysicalNames::VORTICITYZ);
      spSim->addStatsOutputFile(spAvgZ);

      // Create RMS vertical vorticity writer
      IoStats::SharedCartesian1DScalarRMSWriter spRMSZ(new IoStats::Cartesian1DScalarRMSWriter("vorticityz", spAvgZ,  SchemeType::type()));
      spRMSZ->expect(PhysicalNames::VORTICITYZ);
      spSim->addStatsOutputFile(spRMSZ);

     // Create skew vertical vorticity writer
       IoStats::SharedCartesian1DScalarSkewWriter spSkewZ(new IoStats::Cartesian1DScalarSkewWriter("vorticityz", spAvgZ, spRMSZ,  SchemeType::type()));
      spSkewZ->expect(PhysicalNames::VORTICITYZ);
      spSim->addStatsOutputFile(spSkewZ);

      // Create kurt vertical vortcity writer
      IoStats::SharedCartesian1DScalarKurtWriter spKurtZ(new IoStats::Cartesian1DScalarKurtWriter("vorticityz", spAvgZ, spRMSZ,  SchemeType::type()));
      spKurtZ->expect(PhysicalNames::VORTICITYZ);
      spSim->addStatsOutputFile(spKurtZ);

      // streamfunction

      // Create Avg stream function writer
      IoStats::SharedCartesian1DScalarAvgWriter spAvgP(new IoStats::Cartesian1DScalarAvgWriter("streamfunction",SchemeType::type()));
      spAvgP->expect(PhysicalNames::STREAMFUNCTION);
      spSim->addStatsOutputFile(spAvgP);

      // Create RMS stream function writer
      IoStats::SharedCartesian1DScalarRMSWriter spRMSP(new IoStats::Cartesian1DScalarRMSWriter("streamfunction", spAvgP,  SchemeType::type()));
      spRMSP->expect(PhysicalNames::STREAMFUNCTION);
      spSim->addStatsOutputFile(spRMSP);

     // Create skew stream function writer
       IoStats::SharedCartesian1DScalarSkewWriter spSkewP(new IoStats::Cartesian1DScalarSkewWriter("streamfunction", spAvgP, spRMSP,  SchemeType::type()));
      spSkewP->expect(PhysicalNames::STREAMFUNCTION);
      spSim->addStatsOutputFile(spSkewP);

      // Create kurt stream function writer
      IoStats::SharedCartesian1DScalarKurtWriter spKurtP(new IoStats::Cartesian1DScalarKurtWriter("streamfunction", spAvgP, spRMSP,  SchemeType::type()));
      spKurtP->expect(PhysicalNames::STREAMFUNCTION);
      spSim->addStatsOutputFile(spKurtP);

      // velocityx

      // Create Avg velocityx writer
      IoStats::SharedCartesian1DScalarAvgWriter spAvgVx(new IoStats::Cartesian1DScalarAvgWriter("velocityx",SchemeType::type()));
      spAvgVx->expect(PhysicalNames::VELOCITYX);
      spSim->addStatsOutputFile(spAvgVx);

      // velocityy

      // Create Avg velocityy writer
      IoStats::SharedCartesian1DScalarAvgWriter spAvgVy(new IoStats::Cartesian1DScalarAvgWriter("velocityy",SchemeType::type()));
      spAvgVy->expect(PhysicalNames::VELOCITYY);
      spSim->addStatsOutputFile(spAvgVy);

      // thermal dissipation

      // Create Avg thermal dissipation writer
      IoStats::SharedCartesian1DScalarAvgWriter spAvgEth(new IoStats::Cartesian1DScalarAvgWriter("dissTh",SchemeType::type()));
      spAvgEth->expect(PhysicalNames::DISSTH);
      spSim->addStatsOutputFile(spAvgEth);

      // viscous dissipation

      // Create Avg viscous dissipation writer
      IoStats::SharedCartesian1DScalarAvgWriter spAvgEv(new IoStats::Cartesian1DScalarAvgWriter("dissV",SchemeType::type()));
      spAvgEv->expect(PhysicalNames::DISSV);
      spSim->addStatsOutputFile(spAvgEv);
      
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
