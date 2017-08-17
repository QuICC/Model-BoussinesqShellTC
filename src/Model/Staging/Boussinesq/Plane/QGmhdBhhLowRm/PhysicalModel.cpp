/** 
 * @file PhysicalModel.cpp
 * @brief Source of the Boussinesq F-plane with horizontal helicoidal magnetic field QG physical model
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 * @ modified by Stefano Maffei \<maffei.ste@gmail.com\>
 */

/** 
 * Useful information:
 * POLYCOSCOS(a1,k1,a2,k2,a3,k3) = a1*Z^k1 * a2*cos(k2*x) * a3*cos(k3*y)
 * but remember that in the final output x and y are exchanged (transposed)
 * plus, these functions are defined over -1<z<1 , 0<x<2Pi , 0<y<2Pi, despite the actual horizotnal resolution
 * output looks ok in the hdf5 files though (plots with mathematica don't show transposed output)
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
#include MAKE_STR( QUICC_MODEL_PATH/Boussinesq/Plane/QGmhdBhhLowRm/PhysicalModel.hpp )

// Project includes
//
#include MAKE_STR( QUICC_MODEL_PATH/Boussinesq/Plane/QGmhdBhhLowRm/Streamfunction.hpp )
#include MAKE_STR( QUICC_MODEL_PATH/Boussinesq/Plane/QGmhdBhhLowRm/VelocityX.hpp )
#include MAKE_STR( QUICC_MODEL_PATH/Boussinesq/Plane/QGmhdBhhLowRm/VelocityY.hpp )
#include MAKE_STR( QUICC_MODEL_PATH/Boussinesq/Plane/QGmhdBhhLowRm/VelocityZ.hpp )
#include MAKE_STR( QUICC_MODEL_PATH/Boussinesq/Plane/QGmhdBhhLowRm/Transport.hpp )
#include MAKE_STR( QUICC_MODEL_PATH/Boussinesq/Plane/QGmhdBhhLowRm/VorticityZ.hpp )
#include MAKE_STR( QUICC_MODEL_PATH/Boussinesq/Plane/QGmhdBhhLowRm/MeanHeat.hpp )
#include MAKE_STR( QUICC_MODEL_PATH/Boussinesq/Plane/QGmhdBhhLowRm/fbx.hpp )
#include MAKE_STR( QUICC_MODEL_PATH/Boussinesq/Plane/QGmhdBhhLowRm/fby.hpp )
#include MAKE_STR( QUICC_MODEL_PATH/Boussinesq/Plane/QGmhdBhhLowRm/fbz.hpp )
#include MAKE_STR( QUICC_MODEL_PATH/Boussinesq/Plane/QGmhdBhhLowRm/fjz.hpp )
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

namespace QGmhdBhhLowRm {

   const std::string PhysicalModel::PYMODULE = "boussinesq_qgmhdbhhlowrm";

   const std::string PhysicalModel::PYCLASS = "BoussinesqQGmhdBhhLowRm";

   void PhysicalModel::addEquations(SharedSimulation spSim)
   {
      // Add upright streamfunction equation
      spSim->addScalarEquation<Equations::Boussinesq::Plane::QGmhdBhhLowRm::Streamfunction>();

      // Add upright vertical velocity equation
      spSim->addScalarEquation<Equations::Boussinesq::Plane::QGmhdBhhLowRm::VelocityZ>();

      // Add upright transport equation
      spSim->addScalarEquation<Equations::Boussinesq::Plane::QGmhdBhhLowRm::Transport>();

      // Add velocity x
      spSim->addScalarEquation<Equations::Boussinesq::Plane::QGmhdBhhLowRm::VelocityX>();

      // Add velocity y
      spSim->addScalarEquation<Equations::Boussinesq::Plane::QGmhdBhhLowRm::VelocityY>();

      // Add vertical vorticity equation
      spSim->addScalarEquation<Equations::Boussinesq::Plane::QGmhdBhhLowRm::VorticityZ>();

      // Add mean heat computation
      spSim->addScalarEquation<Equations::Boussinesq::Plane::QGmhdBhhLowRm::MeanHeat>();

      // Add fbx computation
      spSim->addScalarEquation<Equations::Boussinesq::Plane::QGmhdBhhLowRm::fbx>();

      // Add fby computation
      spSim->addScalarEquation<Equations::Boussinesq::Plane::QGmhdBhhLowRm::fby>();

      // Add fbz heat computation
      spSim->addScalarEquation<Equations::Boussinesq::Plane::QGmhdBhhLowRm::fbz>();

      // Add fjz computation
      spSim->addScalarEquation<Equations::Boussinesq::Plane::QGmhdBhhLowRm::fjz>();
   }

   void PhysicalModel::addStates(SharedStateGenerator spGen)
   {  
      // ****** Work in progress  ****** //
      // Get flag for the magnetic field to be imposed
      // -1: BX=BY=BZ=0
      // 1: BX=1, BY=BZ=0
      // 2: Helicoidal horizotnal field from Stellmach & Hansen, 2004
//      MHDFloat FB = this->eqParams().nd(NonDimensional::FLAG_BH);
      // ****** ****** //  
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

         // Add BX initial state generation equation
         spExact = spGen->addScalarEquation<Equations::CartesianExactScalarState>();
         spExact->setIdentity(PhysicalNames::BX);
         spExact->setStateType(Equations::CartesianExactStateIds::POLYCOSCOS);
         spExact->setModeOptions(-1e0, 0.0, 1e0, 0.0, 1e0, 0.0);

         // Add BY initial state generation equation
         spExact = spGen->addScalarEquation<Equations::CartesianExactScalarState>();
         spExact->setIdentity(PhysicalNames::BY);
         spExact->setStateType(Equations::CartesianExactStateIds::POLYCOSCOS);
         spExact->setModeOptions(-1e0, 0.0, 1e0, 0.0, 1e0, 0.0);

         // Add BY initial state generation equation
//         spExact = spGen->addScalarEquation<Equations::CartesianExactScalarState>();
//         spExact->setIdentity(PhysicalNames::EMFY);
//         spExact->setStateType(Equations::CartesianExactStateIds::POLYCOSCOS);
//         spExact->setModeOptions(0, 0.0, 0, 0.0, 0, 0.0);
         
         // Add BY initial state generation equation
//         spExact = spGen->addScalarEquation<Equations::CartesianExactScalarState>();
//         spExact->setIdentity(PhysicalNames::EMFX);
//         spExact->setStateType(Equations::CartesianExactStateIds::POLYCOSCOS);
//         spExact->setModeOptions(0, 0.0, 0, 0.0, 0, 0.0);

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
         // An initial condition that satisfies isothermal BCs
//         spExact = spGen->addScalarEquation<Equations::CartesianExactScalarState>();
//         spExact->setIdentity(PhysicalNames::TEMPERATURE);
//         spExact->setStateType(Equations::CartesianExactStateIds::SINCOSCOS);
//         spExact->setModeOptions(1e0, Math::PI, 1e0, 10.0, 1e0, 10.0);

         // Add streamfunction initial state generation equation
         // Random initial condition:
         spRand = spGen->addScalarEquation<Equations::RandomScalarState>();
         spRand->setIdentity(PhysicalNames::STREAMFUNCTION);
         spRand->setSpectrum(-1e-2, 1e-2, 1e4, 1e4, 1e4);
         // ***TEST***
//         spExact = spGen->addScalarEquation<Equations::CartesianExactScalarState>();
//         spExact->setIdentity(PhysicalNames::STREAMFUNCTION);
//         spExact->setStateType(Equations::CartesianExactStateIds::POLYCOSCOS);
//         spExact->setModeOptions(1e0, 2.0, 1e0, 0.159155, 1e0, 0.0);
         // an initial condition that satisfies stress free BCs
//         spExact = spGen->addScalarEquation<Equations::CartesianExactScalarState>();
//         spExact->setIdentity(PhysicalNames::STREAMFUNCTION);
//         spExact->setStateType(Equations::CartesianExactStateIds::COSCOSCOS);
//         spExact->setModeOptions(1e0, Math::PI, 1e0, 10.0, 1e0, 10.0);           

         // Add vertical velocity initial state generation equation
         spRand = spGen->addScalarEquation<Equations::RandomScalarState>();
         spRand->setIdentity(PhysicalNames::VELOCITYZ);
         spRand->setSpectrum(-1e-2, 1e-2, 1e4, 1e4, 1e4);
         // ***TEST***
//         spExact = spGen->addScalarEquation<Equations::CartesianExactScalarState>();
//         spExact->setIdentity(PhysicalNames::VELOCITYZ);
//         spExact->setStateType(Equations::CartesianExactStateIds::POLYCOSCOS);
//         spExact->setModeOptions(1e0, 0.0, 1e0, 0.0, 1e0, 0.0);
         // an intial condition that satisfies non penetration BCs
//         spExact = spGen->addScalarEquation<Equations::CartesianExactScalarState>();
//         spExact->setIdentity(PhysicalNames::VELOCITYZ);
//         spExact->setStateType(Equations::CartesianExactStateIds::SINCOSCOS);
//         spExact->setModeOptions(1e0, Math::PI, 1e0, 10.0, 1e0, 10.0);   

         // Add vertical velocity initial state generation equation
         spExact = spGen->addScalarEquation<Equations::CartesianExactScalarState>();
         spExact->setIdentity(PhysicalNames::DZ_MEANTEMPERATURE);
         spExact->setStateType(Equations::CartesianExactStateIds::POLYCOSCOS);
         spExact->setModeOptions(-1e0, 0.0, 1e0, 0.0, 1e0, 0.0);

         // Add BX initial state generation equation
//         spRand = spGen->addScalarEquation<Equations::RandomScalarState>();
//         spRand->setIdentity(PhysicalNames::BX);
//         spRand->setSpectrum(-1e-2, 1e-2, 1e4, 1e4, 1e4);
         // Helicoidal Bx
         spExact = spGen->addScalarEquation<Equations::CartesianExactScalarState>();
         spExact->setIdentity(PhysicalNames::BX);
         spExact->setStateType(Equations::CartesianExactStateIds::BXHELICOIDAL);
//         spExact->setStateType(Equations::CartesianExactStateIds::CONSTANTFIELD);

         // Add BY initial state generation equation
//         spRand = spGen->addScalarEquation<Equations::RandomScalarState>();
//         spRand->setIdentity(PhysicalNames::BY);
//         spRand->setSpectrum(-1e-2, 1e-2, 1e4, 1e4, 1e4);
         // Helicoidal BY
         spExact = spGen->addScalarEquation<Equations::CartesianExactScalarState>();
         spExact->setIdentity(PhysicalNames::BY);
         spExact->setStateType(Equations::CartesianExactStateIds::BYHELICOIDAL);
//         spExact->setStateType(Equations::CartesianExactStateIds::NULLFIELD);

      }

      // Add output file
      IoVariable::SharedStateFileWriter spOut(new IoVariable::StateFileWriter(SchemeType::type(), SchemeType::isRegular()));
      spOut->expect(PhysicalNames::TEMPERATURE);
      spOut->expect(PhysicalNames::STREAMFUNCTION);
      spOut->expect(PhysicalNames::VELOCITYZ);
      spOut->expect(PhysicalNames::DZ_MEANTEMPERATURE);
      spOut->expect(PhysicalNames::BX);
      spOut->expect(PhysicalNames::BY);
      spGen->addHdf5OutputFile(spOut);
   }

   void PhysicalModel::addVisualizers(SharedVisualizationGenerator spVis)
   {
      // Shared pointer to basic field visualizer
      Equations::SharedScalarFieldVisualizer spField;

      // Add transport field visualization
      spField = spVis->addScalarEquation<Equations::ScalarFieldVisualizer>();
      spField->setFields(true, false);
      spField->setIdentity(PhysicalNames::TEMPERATURE);

      // Add streamfunction field visualization
      spField = spVis->addScalarEquation<Equations::ScalarFieldVisualizer>();
      spField->setFields(true, false);
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
      spField->setFields(true, false);
      spField->setIdentity(PhysicalNames::VELOCITYZ);

      // Add background temperature profile visualization
      spField = spVis->addScalarEquation<Equations::ScalarFieldVisualizer>();
      spField->setFields(true, false);
      spField->setIdentity(PhysicalNames::DZ_MEANTEMPERATURE);

      // Add background temperature profile visualization
      spField = spVis->addScalarEquation<Equations::ScalarFieldVisualizer>();
      spField->setFields(true, false);
      spField->setIdentity(PhysicalNames::BX);

      // Add background temperature profile visualization
      spField = spVis->addScalarEquation<Equations::ScalarFieldVisualizer>();
      spField->setFields(true, false);
      spField->setIdentity(PhysicalNames::BY);

      // Add background temperature profile visualization
      spField = spVis->addScalarEquation<Equations::ScalarFieldVisualizer>();
      spField->setFields(true, false);
      spField->setIdentity(PhysicalNames::VORTICITYZ);

      // Add background temperature profile visualization
      spField = spVis->addScalarEquation<Equations::ScalarFieldVisualizer>();
      spField->setFields(true, false);
      spField->setIdentity(PhysicalNames::FBX);

      // Add background temperature profile visualization
      spField = spVis->addScalarEquation<Equations::ScalarFieldVisualizer>();
      spField->setFields(true, false);
      spField->setIdentity(PhysicalNames::FBY);

      // Add background temperature profile visualization
      spField = spVis->addScalarEquation<Equations::ScalarFieldVisualizer>();
      spField->setFields(true, false);
      spField->setIdentity(PhysicalNames::FBZ);

      // Add background temperature profile visualization
      spField = spVis->addScalarEquation<Equations::ScalarFieldVisualizer>();
      spField->setFields(true, false);
      spField->setIdentity(PhysicalNames::FJZ);

      // Add output file
      IoVariable::SharedVisualizationFileWriter spOut(new IoVariable::VisualizationFileWriter(SchemeType::type()));
      spOut->expect(PhysicalNames::TEMPERATURE);
      spOut->expect(PhysicalNames::STREAMFUNCTION);
      spOut->expect(PhysicalNames::VELOCITYX);
      spOut->expect(PhysicalNames::VELOCITYY);
      spOut->expect(PhysicalNames::VELOCITYZ);
      spOut->expect(PhysicalNames::DZ_MEANTEMPERATURE);
      spOut->expect(PhysicalNames::BX);
      spOut->expect(PhysicalNames::BY);
      spOut->expect(PhysicalNames::VORTICITYZ);
      spOut->expect(PhysicalNames::FBX);
      spOut->expect(PhysicalNames::FBY);
      spOut->expect(PhysicalNames::FBZ);
      spOut->expect(PhysicalNames::FJZ);
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
      spIn->expect(PhysicalNames::BX);
      spIn->expect(PhysicalNames::BY);
      spIn->expect(PhysicalNames::VORTICITYZ);
      spIn->expect(PhysicalNames::FBY);
      spIn->expect(PhysicalNames::FBX);
      spIn->expect(PhysicalNames::FBZ);
      spIn->expect(PhysicalNames::FJZ);

      // Set simulation state
      spVis->setInitialState(spIn);
   }

   void PhysicalModel::addAsciiOutputFiles(SharedSimulation spSim)
   {
      // Create Nusselt number writer
      //IoVariable::SharedNusseltWriter spNusselt(new IoVariable::NusseltWriter(SchemeType::type()));
      //spNusselt->expect(PhysicalNames::DZ_MEANTEMPERATURE);
      //spSim->addAsciiOutputFile(spNusselt);

      // Create Nusselt number writer
      IoVariable::SharedCartesian1DNusseltZWriter spNusselt(new IoVariable::Cartesian1DNusseltZWriter(SchemeType::type()));
      spNusselt->expect(PhysicalNames::DZ_MEANTEMPERATURE);
      spSim->addAsciiOutputFile(spNusselt);

      // Create temperature energy writer
      // NB: It's energy per unit volume (is it?)
      IoVariable::SharedCartesian1DScalarEnergyWriter spTemp(new IoVariable::Cartesian1DScalarEnergyWriter("temperature", SchemeType::type()));
      spTemp->expect(PhysicalNames::TEMPERATURE);
      spSim->addAsciiOutputFile(spTemp);

      // Create kinetic energy writer
      // NB: It's kinetic energy per unit volume
      IoVariable::SharedCartesian1DStreamEnergyWriter spStream(new IoVariable::Cartesian1DStreamEnergyWriter("kinetic", SchemeType::type()));
      spStream->expect(PhysicalNames::STREAMFUNCTION);
      spStream->expect(PhysicalNames::VELOCITYZ);
      spSim->addAsciiOutputFile(spStream);

      // Create magnetic energy writer
      // NB: It's mangetic energy per unit volume
      IoVariable::SharedCartesian1DMagneticEnergyWriter spMag(new IoVariable::Cartesian1DMagneticEnergyWriter("magnetic", SchemeType::type()));
      spMag->expect(PhysicalNames::BX);
      spMag->expect(PhysicalNames::BY);
      spSim->addAsciiOutputFile(spMag);

      // Create fluctuating magnetic energy writer
      // NB: It's mangetic energy per unit volume
      IoVariable::SharedCartesian1DFluctuatingMagneticEnergyWriter spFlMag(new IoVariable::Cartesian1DFluctuatingMagneticEnergyWriter("fluct_magnetic", SchemeType::type()));
      spFlMag->expect(PhysicalNames::FBX);
      spFlMag->expect(PhysicalNames::FBY);
      spFlMag->expect(PhysicalNames::FBZ);
      spSim->addAsciiOutputFile(spFlMag);

      // Create cartesian kinetic energy writer
      // to be compared with kinetic_energy, created by calculating velocityx and velocityy
      // NB: It's mangetic energy per unit volume
      IoVariable::SharedCartesian1DKineticCartesianWriter spCartK(new IoVariable::Cartesian1DKineticCartesianWriter("cartesian_kinetic", SchemeType::type()));
      spCartK->expect(PhysicalNames::VELOCITYX);
      spCartK->expect(PhysicalNames::VELOCITYY);
      spCartK->expect(PhysicalNames::VELOCITYZ);
      spSim->addAsciiOutputFile(spCartK);

      // Create temperature energy writer
      //   IoVariable::SharedCartesian1DScalarEnergyWriter spMag(new IoVariable::Cartesian1DScalarEnergyWriter("magX", SchemeType::type()));
      //   spMag->expect(PhysicalNames::BX);
      //   spSim->addAsciiOutputFile(spMag);

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

      if(false)
      {

      // velocityx

      // Create Avg velocityx writer
      IoStats::SharedCartesian1DScalarAvgWriter spAvgVx(new IoStats::Cartesian1DScalarAvgWriter("velocityx",SchemeType::type()));
      spAvgVx->expect(PhysicalNames::VELOCITYX);
      spSim->addStatsOutputFile(spAvgVx);

      // Create RMS velocityx  writer
      IoStats::SharedCartesian1DScalarRMSWriter spRMSVx(new IoStats::Cartesian1DScalarRMSWriter("velocityx", spAvgVx,  SchemeType::type()));
      spRMSVx->expect(PhysicalNames::VELOCITYX);
      spSim->addStatsOutputFile(spRMSVx);

      // Create skew velocityx  writer       
      IoStats::SharedCartesian1DScalarSkewWriter spSkewVx(new IoStats::Cartesian1DScalarSkewWriter("velocityx", spAvgVx, spRMSVx,  SchemeType::type()));
      spSkewVx->expect(PhysicalNames::VELOCITYX);
      spSim->addStatsOutputFile(spSkewVx);

      // Create kurt velocityx  writer
      IoStats::SharedCartesian1DScalarKurtWriter spKurtVx(new IoStats::Cartesian1DScalarKurtWriter("velocityx", spAvgVx, spRMSVx,  SchemeType::type()));
      spKurtVx->expect(PhysicalNames::VELOCITYX);
      spSim->addStatsOutputFile(spKurtVx);

      // velocityy

      // Create Avg velocityy writer
      IoStats::SharedCartesian1DScalarAvgWriter spAvgVy(new IoStats::Cartesian1DScalarAvgWriter("velocityy",SchemeType::type()));
      spAvgVy->expect(PhysicalNames::VELOCITYY);
      spSim->addStatsOutputFile(spAvgVy);

      // Create RMS velocityy  writer
      IoStats::SharedCartesian1DScalarRMSWriter spRMSVy(new IoStats::Cartesian1DScalarRMSWriter("velocityy", spAvgVy,  SchemeType::type()));
      spRMSVy->expect(PhysicalNames::VELOCITYY);
      spSim->addStatsOutputFile(spRMSVy);

     // Create skew velocityy  writer       
     IoStats::SharedCartesian1DScalarSkewWriter spSkewVy(new IoStats::Cartesian1DScalarSkewWriter("velocityy", spAvgVy, spRMSVy,  SchemeType::type()));
      spSkewVy->expect(PhysicalNames::VELOCITYY);
      spSim->addStatsOutputFile(spSkewVy);

      // Create kurt velocityy  writer
      IoStats::SharedCartesian1DScalarKurtWriter spKurtVy(new IoStats::Cartesian1DScalarKurtWriter("velocityy", spAvgVy, spRMSVy,  SchemeType::type()));
      spKurtVy->expect(PhysicalNames::VELOCITYY);
      spSim->addStatsOutputFile(spKurtVy);

      // fjz

      // Create Avg vertical current writer
      IoStats::SharedCartesian1DScalarAvgWriter spAvgFJZ(new IoStats::Cartesian1DScalarAvgWriter("fjz",SchemeType::type()));
      spAvgFJZ->expect(PhysicalNames::FJZ);
      spSim->addStatsOutputFile(spAvgFJZ);

      // Create RMS vertical current writer
      IoStats::SharedCartesian1DScalarRMSWriter spRMSFJZ(new IoStats::Cartesian1DScalarRMSWriter("fjz", spAvgFJZ,  SchemeType::type()));
      spRMSFJZ->expect(PhysicalNames::FJZ);
      spSim->addStatsOutputFile(spRMSFJZ);

     // Create skew vertical current writer
       IoStats::SharedCartesian1DScalarSkewWriter spSkewFJZ(new IoStats::Cartesian1DScalarSkewWriter("fjz", spAvgFJZ, spRMSFJZ,  SchemeType::type()));
      spSkewFJZ->expect(PhysicalNames::FJZ);
      spSim->addStatsOutputFile(spSkewFJZ);

      // Create kurt vertical current writer
      IoStats::SharedCartesian1DScalarKurtWriter spKurtFJZ(new IoStats::Cartesian1DScalarKurtWriter("fjz", spAvgFJZ, spRMSFJZ,  SchemeType::type()));
      spKurtFJZ->expect(PhysicalNames::FJZ);
      spSim->addStatsOutputFile(spKurtFJZ);

      }

      // fbz

      // Create Avg fluctuating vertical magnetic field writer
      IoStats::SharedCartesian1DScalarAvgWriter spAvgFBZ(new IoStats::Cartesian1DScalarAvgWriter("fbz",SchemeType::type()));
      spAvgFBZ->expect(PhysicalNames::FBZ);
      spSim->addStatsOutputFile(spAvgFBZ);

      // Create RMS fluctuating vertical magnetic field  writer
      IoStats::SharedCartesian1DScalarRMSWriter spRMSFBZ(new IoStats::Cartesian1DScalarRMSWriter("fbz", spAvgFBZ,  SchemeType::type()));
      spRMSFBZ->expect(PhysicalNames::FBZ);
      spSim->addStatsOutputFile(spRMSFBZ);

     // Create skew fluctuating vertical magnetic field  writer
       IoStats::SharedCartesian1DScalarSkewWriter spSkewFBZ(new IoStats::Cartesian1DScalarSkewWriter("fbz", spAvgFBZ, spRMSFBZ,  SchemeType::type()));
      spSkewFBZ->expect(PhysicalNames::FBZ);
      spSim->addStatsOutputFile(spSkewFBZ);

      // Create kurt  fluctuating vertical magnetic field  writer
      IoStats::SharedCartesian1DScalarKurtWriter spKurtFBZ(new IoStats::Cartesian1DScalarKurtWriter("fbz", spAvgFBZ, spRMSFBZ,  SchemeType::type()));
      spKurtFBZ->expect(PhysicalNames::FBZ);
      spSim->addStatsOutputFile(spKurtFBZ);

      // fbx

      // Create Avg fluctuating x magnetic field writer
      IoStats::SharedCartesian1DScalarAvgWriter spAvgFBX(new IoStats::Cartesian1DScalarAvgWriter("fbx",SchemeType::type()));
      spAvgFBX->expect(PhysicalNames::FBX);
      spSim->addStatsOutputFile(spAvgFBX);

      // Create RMS fluctuating x magnetic field  writer
      IoStats::SharedCartesian1DScalarRMSWriter spRMSFBX(new IoStats::Cartesian1DScalarRMSWriter("fbx", spAvgFBX,  SchemeType::type()));
      spRMSFBX->expect(PhysicalNames::FBX);
      spSim->addStatsOutputFile(spRMSFBX);

      // Create skew fluctuating x magnetic field  writer
      IoStats::SharedCartesian1DScalarSkewWriter spSkewFBX(new IoStats::Cartesian1DScalarSkewWriter("fbx", spAvgFBX, spRMSFBX,  SchemeType::type()));
      spSkewFBX->expect(PhysicalNames::FBX);
      spSim->addStatsOutputFile(spSkewFBX);

      // Create kurt fluctuating x magnetic field  writer
      IoStats::SharedCartesian1DScalarKurtWriter spKurtFBX(new IoStats::Cartesian1DScalarKurtWriter("fbx", spAvgFBX, spRMSFBX,  SchemeType::type()));
      spKurtFBX->expect(PhysicalNames::FBX);
      spSim->addStatsOutputFile(spKurtFBX);

      // fby

      // Create Avg fluctuating y magnetic field writer
      IoStats::SharedCartesian1DScalarAvgWriter spAvgFBY(new IoStats::Cartesian1DScalarAvgWriter("fby",SchemeType::type()));
      spAvgFBY->expect(PhysicalNames::FBY);
      spSim->addStatsOutputFile(spAvgFBY);

      // Create RMS fluctuating y magnetic field  writer
      IoStats::SharedCartesian1DScalarRMSWriter spRMSFBY(new IoStats::Cartesian1DScalarRMSWriter("fby", spAvgFBY,  SchemeType::type()));
      spRMSFBY->expect(PhysicalNames::FBY);
      spSim->addStatsOutputFile(spRMSFBY);

     // Create skew fluctuating y magnetic field  writer       
      IoStats::SharedCartesian1DScalarSkewWriter spSkewFBY(new IoStats::Cartesian1DScalarSkewWriter("fby", spAvgFBY, spRMSFBY,  SchemeType::type()));
      spSkewFBY->expect(PhysicalNames::FBY);
      spSim->addStatsOutputFile(spSkewFBY);

      // Create kurt fluctuating y magnetic field  writer
      IoStats::SharedCartesian1DScalarKurtWriter spKurtFBY(new IoStats::Cartesian1DScalarKurtWriter("fby", spAvgFBY, spRMSFBY,  SchemeType::type()));
      spKurtFBY->expect(PhysicalNames::FBY);
      spSim->addStatsOutputFile(spKurtFBY);


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
      spState->expect(PhysicalNames::BY);
      spState->expect(PhysicalNames::BX);
      spState->expect(PhysicalNames::FBX);
      spState->expect(PhysicalNames::FBY);
      spState->expect(PhysicalNames::FBZ);
      spState->expect(PhysicalNames::FJZ);


      spSim->addHdf5OutputFile(spState);
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

      spInit->expect(PhysicalNames::BY);
      spInit->expect(PhysicalNames::BX);
//      spInit->expect(PhysicalNames::PRESSURE);

      // Set simulation state
      spSim->setInitialState(spInit);
   }
}
}
}
}
}
