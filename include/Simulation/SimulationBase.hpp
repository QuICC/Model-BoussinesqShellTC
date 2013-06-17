/** \file Simulation.hpp
 *  \brief High level implementation of a base for the simulations
 */

#ifndef SIMULATIONBASE_HPP
#define SIMULATIONBASE_HPP

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
#include "SparseSolvers/SparseTrivialCoordinator.hpp"
#include "SparseSolvers/SparseLinearCoordinator.hpp"
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
    * @brief High level implementation of a base for the simulations
    */
   class SimulationBase
   {
      public:
         /// Typedef for a shared scalar equation iterator
         typedef std::vector<Equations::SharedIScalarEquation>::iterator   ScalarEquation_iterator;

         /// Typedef for a shared vector equation iterator
         typedef std::vector<Equations::SharedIVectorEquation>::iterator   VectorEquation_iterator;

         /// Typedef for a shared scalar equation range
         typedef std::pair<ScalarEquation_iterator, ScalarEquation_iterator>  ScalarEquation_range;

         /// Typedef for a shared vector equation range
         typedef std::pair<VectorEquation_iterator, VectorEquation_iterator>  VectorEquation_range;

         /**
          * @brief Constructor
          */
         SimulationBase();

         /**
          * @brief Simple empty destructor
          */
         virtual ~SimulationBase();

         /**
          * @brief Initialise the base components of the simulation
          */
         void initBase();

         /**
          * @brief Initialise the different components of the simulation
          *
          * @param spBcs Boundary condition information
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
         template <typename TEquation> SharedPtrMacro<TEquation> addScalarEquation();

         /**
          * @brief Add vector equation to solver
          */
         template <typename TEquation> SharedPtrMacro<TEquation> addVectorEquation();

         /**
          * @brief Set the base simulation configuration file and parameters
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
         /**
          * @brief Initialise the solvers (done just before preRun)
          */
         virtual void initSolvers();

         /**
          * @brief Compute the nonlinear terms
          */
         void computeNonlinear();

         /**
          * @brief Solve the trivial equations
          */
         void solveTrivialEquations();

         /**
          * @brief Solve the diagnostic equations
          */
         void solveDiagnosticEquations();

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
          * @brief Storage for the range of scalar prognostic equations
          */
         ScalarEquation_range mScalarPrognosticRange;

         /**
          * @brief Storage for the range of vector prognostic equations
          */
         VectorEquation_range mVectorPrognosticRange;

         /**
          * @brief Storage for the range of scalar solver equations
          */
         ScalarEquation_range mScalarDiagnosticRange;

         /**
          * @brief Storage for the range of vector solver equations
          */
         VectorEquation_range mVectorDiagnosticRange;

         /**
          * @brief Storage for the range of scalar trivial equations
          */
         ScalarEquation_range mScalarTrivialRange;

         /**
          * @brief Storage for the range of vectort trivial equations
          */
         VectorEquation_range mVectorTrivialRange;

         /**
          * @brief Simulation run control
          */
         SimulationRunControl mSimRunCtrl;

         /**
          * @brief Trivial solver coordinator
          */
         Solver::SparseTrivialCoordinator mTrivialCoordinator;

         /**
          * @brief Linear solver coordinator
          */
         Solver::SparseLinearCoordinator mLinearCoordinator;

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

      private:
         /**
          * @brief Add addition configuration file parts
          *
          * @param spCfgFile  Configuration file
          */
         virtual void addConfigurationPart(IoConfig::SharedConfigurationReader spCfgFile);

         /**
          * @brief Initialise the implementation specific base components
          */
         virtual void initAdditionalBase();

         /**
          * @brief Do operations required just before starting the main loop
          */
         virtual void preRun() = 0;

         /**
          * @brief Do operations required during the main loop
          */
         virtual void mainRun() = 0;

         /**
          * @brief Do operations required just after finishing the main loop
          */
         virtual void postRun() = 0;

         /**
          * @brief Initialise the variables required by the simulation
          *
          * @param varInfo Global variable requirements
          */
         void initVariables(VariableRequirement& varInfo);

         /**
          * @brief Map the variables to the equations
          *
          * @param nonInfo Global nonlinear requirements
          */
         void mapEquationVariables(std::set<PhysicalNames::Id>& nonInfo);

         /**
          * @brief Initialise the transform coordinator
          *
          * @param varInfo Global variable requirements
          * @param nonInfo Global nonlinear requirements
          */
         void initTransformCoordinator(const VariableRequirement& varInfo, const std::set<PhysicalNames::Id>& nonInfo);

         /**
          * @brief Initialise the equations (generate operators, etc)
          *
          * @param spBcs Boundary condition information
          */
         void setupEquations(const SharedSimulationBoundary spBcs);

         /**
          * @brief Sort equations and store information for timestep/solver/nothing ranges
          */
         void sortEquations();

         /**
          * @brief Setup the output files added by the model
          */
         void setupOutput();

         /**
          * @brief Allow for implementation specific output tuning
          */
         virtual void tuneOutput();
   };

   /**
    * @brief Compute the scalar equation type flag for time/solver/trivial ordering
    */
   int computeScalarEquationType(Equations::SharedIScalarEquation eqA);

   /**
    * @brief Compute the vector equation type flag for time/solver/trivial ordering
    */
   int computeVectorEquationType(Equations::SharedIVectorEquation eqA);

   /**
    * @brief Sorting function for the scalar equations to obtain time/solver/trivial ordering
    */
   bool sortScalarEquationType(Equations::SharedIScalarEquation eqA, Equations::SharedIScalarEquation eqB);

   /**
    * @brief Sorting function for the vector equations to obtain time/solver/trivial ordering
    */
   bool sortVectorEquationType(Equations::SharedIVectorEquation eqA, Equations::SharedIVectorEquation eqB);

   template <typename TScheme> void SimulationBase::initResolution()
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

   template <typename TModel> SharedSimulationBoundary SimulationBase::createBoundary()
   {
      return TModel::createBoundary(this->mSimIoCtrl.configBoundary());
   }

   template <int DIMENSION> void SimulationBase::setConfiguration(const std::string& type, const std::vector<bool>& isPeriodicBox, const std::vector<std::string>& bcNames, const std::vector<std::string>& ndNames)
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

      // Add addition configuration parts
      this->addConfigurationPart(spCfgFile);

      // Set configuration file
      this->mSimIoCtrl.setConfigurationFile(spCfgFile);
   }

   template <typename TEquation> SharedPtrMacro<TEquation> SimulationBase::addScalarEquation()
   {
      // Create shared scalar equation
      SharedPtrMacro<TEquation>  spEq(new TEquation(this->mspEqParams));

      // Add share scalar equation
      this->mScalarEquations.push_back(spEq);

      return spEq;
   }

   template <typename TEquation> SharedPtrMacro<TEquation> SimulationBase::addVectorEquation()
   {
      // Create shared scalar equation
      SharedPtrMacro<TEquation>  spEq(new TEquation(this->mspEqParams));

      // Add share scalar equation
      this->mVectorEquations.push_back(spEq);

      return spEq;
   }

   /// Typedef for a shared pointer of a SimulationBase
   typedef SharedPtrMacro<SimulationBase> SharedSimulationBase;

}

#endif // SIMULATIONBASE_HPP
