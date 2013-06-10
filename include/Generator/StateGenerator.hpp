/** \file StateGenerator.hpp
 *  \brief High level implementation of a general state generator
 */

#ifndef STATEGENERATOR_HPP
#define STATEGENERATOR_HPP

// Configuration includes
//
#include "SmartPointers/SharedPtrMacro.h"

// System includes
//

// External includes
//

// Project includes
//
#include "Simulation/SimulationIoControl.hpp"
#include "Equations/EquationParameters.hpp"
#include "Equations/IScalarDEquation.hpp"
#include "Equations/IVectorDEquation.hpp"
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

namespace GeoMHDiSCC {

   /**
    * @brief High level implementation of a general state generator
    */
   class StateGenerator
   {
      public:
         /**
          * @brief Constructor
          */
         StateGenerator();

         /**
          * @brief Simple empty destructor
          */
         ~StateGenerator();

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
          * @brief Add scalar state equation
          */
         template <typename TEquation> void addScalarState(const PhysicalNames::Id name);

         /**
          * @brief Add vector state equation
          */
         template <typename TEquation> void addVectorState();

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

      protected:

      private:
         /**
          * @brief Write the output if required
          */
         void writeOutput();

         /**
          * @brief Initialise the variables required by the simulation
          */
         void initVariables();

         /**
          * @brief Initialise and setup the equations added by the model
          */
         void setupEquations();

         /**
          * @brief Setup the output files adde by the model
          */
         void setupOutput();

         /**
          * @brief Shared resolution
          */
         SharedResolution mspRes;

         /**
          * @brief Shared resolution
          */
         Equations::SharedEquationParameters mspEqParams;

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
         std::vector<Equations::SharedIScalarDEquation> mScalarEquations;

         /**
          * @brief Storage for vector equations
          */
         std::vector<Equations::SharedIVectorDEquation> mVectorEquations;

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
   };

   template <typename TScheme> void StateGenerator::initResolution()
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

   template <typename TModel> SharedSimulationBoundary StateGenerator::createBoundary()
   {
      return TModel::createBoundary(this->mSimIoCtrl.configBoundary());
   }

   template <int DIMENSION> void StateGenerator::setConfiguration(const std::string& type, const std::vector<bool>& isPeriodicBox, const std::vector<std::string>& bcNames, const std::vector<std::string>& ndNames)
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

   template <typename TEquation> void StateGenerator::addScalarState(const PhysicalNames::Id name)
   {
      // Create shared scalar equation
      SharedPtrMacro<TEquation>  spEq(new TEquation(this->mspEqParams, name));

      // Add share scalar equation
      this->mScalarEquations.push_back(spEq);
   }

   template <typename TEquation> void StateGenerator::addVectorState()
   {
      // Create shared scalar equation
      SharedPtrMacro<TEquation>  spEq(new TEquation(this->mspEqParams));

      // Add share scalar equation
      this->mVectorEquations.push_back(spEq);
   }

   /// Typedef for a shared pointer of a StateGenerator
   typedef SharedPtrMacro<StateGenerator> SharedStateGenerator;

}

#endif // STATEGENERATOR_HPP
