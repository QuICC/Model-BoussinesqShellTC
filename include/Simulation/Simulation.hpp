/** \file Simulation.hpp
 *  \brief High level implementation of a simulation
 */

#ifndef SIMULATION_HPP
#define SIMULATION_HPP

// Configuration includes
//
#include "SmartPointers/SharedPtrMacro.h"

// System includes
//

// External includes
//

// Project includes
//
#include "Timers/ExecutionTimer.hpp"
#include "Simulation/SimulationRunControl.hpp"
#include "Simulation/SimulationIoControl.hpp"
#include "Simulation/SimulationBoundary.hpp"
#include "Equations/EquationParameters.hpp"
#include "Equations/IScalarEquation.hpp"
#include "Equations/IVectorEquation.hpp"
#include "Timesteppers/Timestepper.hpp"
#include "IoConfig/ConfigurationReader.hpp"
#include "TypeSelectors/TransformSelector.hpp"
#include "TypeSelectors/VariableSelector.hpp"
#include "TypeSelectors/ParallelSelector.hpp"
#include "TransformGroupers/IForwardGrouper.hpp"
#include "TransformGroupers/IBackwardGrouper.hpp"
#include "LoadSplitter/LoadSplitter.hpp"
#include "IoConfig/ConfigParts/PhysicalPart.hpp"
#include "IoConfig/ConfigParts/BoundaryPart.hpp"
#include "IoVariable/IVariableHdf5Reader.hpp"
#include "Diagnostics/DiagnosticCoordinator.hpp"

namespace GeoMHDiSCC {

   /**
    * @brief High level implementation of a simulation's execution steps.
    */
   class Simulation
   {
      public:
         /// Typedef for a shared scalar equation iterator
         typedef std::vector<Equations::SharedIScalarEquation>::iterator   ScalarEquationIteratorType;

         /// Typedef for a shared vector equation iterator
         typedef std::vector<Equations::SharedIVectorEquation>::iterator   VectorEquationIteratorType;

         /// Typedef for a shared scalar equation range
         typedef std::pair<ScalarEquationIteratorType, ScalarEquationIteratorType>  ScalarEquationRangeType;

         /// Typedef for a shared vector equation range
         typedef std::pair<VectorEquationIteratorType, VectorEquationIteratorType>  VectorEquationRangeType;

         /**
          * @brief Constructor
          */
         Simulation();

         /**
          * @brief Simple empty destructor
          */
         ~Simulation();

         /**
          * @brief Initialise the base components of the simulation
          */
         void initBase();

         /**
          * @brief Initialise the different components of the simulation
          */
         void init(const SharedSimulationBoundary spBcs);

         /**
          * @brief Run the simulation
          */
         void run();

         /**
          * @brief Finalise simulation run
          */
         void finalize();

         /**
          * @brief Initialise the resolution
          */
         template <typename TScheme> void initResolution();

         /**
          * @brief Create the simulation wide boundary conditions
          */
         template <typename TModel> SharedSimulationBoundary createBoundary();

         /**
          * @brief Add scalar equation to solver
          */
         template <typename TEquation> void addScalarEquation();

         /**
          * @brief Add vector equation to solver
          */
         template <typename TEquation> void addVectorEquation();

         /**
          * @brief Set the simulation configuration file and parameters
          *
          * @param bcNames Vector of names for the boundary conditions
          */
         template <int DIMENSION> void setConfiguration(const std::string& type, const std::vector<bool>& isPeriodicBox, const std::vector<std::string>& bcNames, const std::vector<std::string>& ndNames);

         /**
          * @brief Set initial state through input file
          *
          * @param spInitFile Shared initial state file
          */
         void setInitialState(IoVariable::SharedIVariableHdf5Reader spInitFile);

         /**
          * @brief Add ASCII output file to solver
          *
          * @param spOutFile Shared ASCII output file
          *
          * \mhdBug Fake declaration
          */
         void addOutputFile(int spOutFile);//SharedAscii spOutFile);

         /**
          * @brief Add HDF5 output file to solver
          *
          * @param spOutFile Shared HDF5 output file
          */
         void addOutputFile(IoVariable::SharedIVariableHdf5NWriter spOutFile);

      protected:

      private:
         /**
          * @brief Do operations required just before starting the time integration
          */
         void preRun();

         /**
          * @brief Compute the nonlinear terms
          */
         void computeNonlinear();

         /**
          * @brief Timestep the equations
          */
         void timestepEquations();

         /**
          * @brief Write the output if required
          */
         void writeOutput();

         /**
          * @brief Do operations required just after finishing the time integration
          */
         void postRun();

         /**
          * @brief Initialise the variables required by the simulation
          */
         void initVariables();

         /**
          * @brief Map the variables to the equations
          */
         void mapEquationVariables();

         /**
          * @brief Initialise the equations (generate operators, etc)
          */
         void setupEquations(const SharedSimulationBoundary spBcs);

         /**
          * @brief Sort equations and store information for timestep/solver/nothing ranges
          */
         void sortEquations();

         /**
          * @brief Setup the output files adde by the model
          */
         void setupOutput();

         /**
          * @brief Execution timer
          */
         ExecutionTimer mExecutionTimer;

         /**
          * @brief Shared resolution
          */
         SharedResolution mspRes;

         /**
          * @brief Shared resolution
          */
         Equations::SharedEquationParameters mspEqParams;

         /**
          * @brief Simulation run control
          */
         SimulationRunControl mSimRunCtrl;

         /**
          * @brief Timestepper
          */
         Timestep::Timestepper mTimestepper;

         /**
          * @brief Simulation IO control
          */
         SimulationIoControl mSimIoCtrl;

         /**
          * @brief Transform coordinator
          */
         Transform::TransformCoordinatorType mTransformCoordinator;

         /**
          * @brief Storage for scalar equations
          */
         std::vector<Equations::SharedIScalarEquation> mScalarEquations;

         /**
          * @brief Storage for vector equations
          */
         std::vector<Equations::SharedIVectorEquation> mVectorEquations;

         /**
          * @brief Storage for the range of scalar timestepping equations
          */
         ScalarEquationRangeType mScalarPrognosticRange;

         /**
          * @brief Storage for the range of vector timestepping equations
          */
         VectorEquationRangeType mVectorPrognosticRange;

         /**
          * @brief Storage for the range of scalar solver equations
          */
         ScalarEquationRangeType mScalarDiagnosticRange;

         /**
          * @brief Storage for the range of vector solver equations
          */
         VectorEquationRangeType mVectorDiagnosticRange;

         /**
          * @brief Storage for the range of scalar "direct" equations
          */
         ScalarEquationRangeType mScalarDirectRange;

         /**
          * @brief Storage for the range of vector "direct" equations
          */
         VectorEquationRangeType mVectorDirectRange;

         /**
          * @brief Map between name and pointer for the scalar variables
          */
         std::map<PhysicalNames::Id, Datatypes::SharedScalarVariableType>  mScalarVariables;

         /**
          * @brief Map between name and pointer for the vector variables
          */
         std::map<PhysicalNames::Id, Datatypes::SharedVectorVariableType>  mVectorVariables;

         /**
          * @brief Storage for a shared forward transform grouper
          */
         Transform::SharedIForwardGrouper   mspFwdGrouper;

         /**
          * @brief Storage for a shared backward transform grouper
          */
         Transform::SharedIBackwardGrouper   mspBwdGrouper;

         /**
          * @brief Diagnostic coordinator
          */
         Diagnostics::DiagnosticCoordinator  mDiagnostics;
   };

   /**
    * @brief Compute the scalar equation type flag for time/solver/directo ordering
    */
   int computeScalarEquationType(Equations::SharedIScalarEquation eqA);

   /**
    * @brief Compute the vector equation type flag for time/solver/directo ordering
    */
   int computeVectorEquationType(Equations::SharedIVectorEquation eqA);

   /**
    * @brief Sorting function for the scalar equations to obtain time/solver/directo ordering
    */
   bool sortScalarEquationType(Equations::SharedIScalarEquation eqA, Equations::SharedIScalarEquation eqB);

   /**
    * @brief Sorting function for the vector equations to obtain time/solver/directo ordering
    */
   bool sortVectorEquationType(Equations::SharedIVectorEquation eqA, Equations::SharedIVectorEquation eqB);

   template <typename TScheme> void Simulation::initResolution()
   {
      // Create the load splitter
      Parallel::LoadSplitter splitter(FrameworkMacro::id(), FrameworkMacro::nCpu());

      // Extract dimensions from configuration file
      ArrayI dim = this->mSimIoCtrl.configDimension();

      // Initialise the load splitter
      splitter.init<TScheme>(dim);

      // Get best splitting resolution object
      std::pair<SharedResolution, Parallel::SplittingDescription>  best = splitter.bestSplitting();

      // Store the shared resolution object
      this->mspRes = best.first;

      // Extract box scale from configuration file
      Array box = this->mSimIoCtrl.configBoxScale();

      // Set the box scale
      this->mspRes->setBoxScale(box);

      // Initialise the transform grouper
      Parallel::setGrouper(best.second, this->mspFwdGrouper, this->mspBwdGrouper);
   }

   template <typename TModel> SharedSimulationBoundary Simulation::createBoundary()
   {
      return TModel::createBoundary(this->mSimIoCtrl.configBoundary());
   }

   template <int DIMENSION> void Simulation::setConfiguration(const std::string& type, const std::vector<bool>& isPeriodicBox, const std::vector<std::string>& bcNames, const std::vector<std::string>& ndNames)
   {
      // Create shared configuration file
      IoConfig::SharedConfigurationReader spCfgFile(new IoConfig::ConfigurationReader(DIMENSION, isPeriodicBox, type));

      // Create the equation parameter shared pointer
      Equations::SharedEquationParameters spEqParams(new Equations::EquationParameters);
      this->mspEqParams = spEqParams;

      // Create the equation parameter dependent configuration part
      IoConfig::SharedPhysicalPart   spPhys(new IoConfig::PhysicalPart(ndNames));

      // Add physical part to configuration file
      spCfgFile->addPart(IoConfig::SimulationBlocks::PHYSICAL, spPhys);

      // Create the boundary condition configuration part
      IoConfig::SharedBoundaryPart   spBound(new IoConfig::BoundaryPart(bcNames));

      // Add boundary part to configuration file
      spCfgFile->addPart(IoConfig::SimulationBlocks::BOUNDARY, spBound);

      // Set configuration file
      this->mSimIoCtrl.setConfigurationFile(spCfgFile);
   }

   template <typename TEquation> void Simulation::addScalarEquation()
   {
      // Create shared scalar equation
      SharedPtrMacro<TEquation>  spEq(new TEquation(this->mspEqParams));

      // Add share scalar equation
      this->mScalarEquations.push_back(spEq);
   }

   template <typename TEquation> void Simulation::addVectorEquation()
   {
      // Create shared scalar equation
      SharedPtrMacro<TEquation>  spEq(new TEquation(this->mspEqParams));

      // Add share scalar equation
      this->mVectorEquations.push_back(spEq);
   }

   /// Typedef for a shared pointer of a Simulation
   typedef SharedPtrMacro<Simulation> SharedSimulation;

}

#endif // SIMULATION_HPP
