/** \file Simulation.cpp
 *  \brief Source of the high level simulation
 */

// Configuration includes
//
#include "Profiler/ProfilerMacro.h"
#include "StorageProfiler/StorageProfilerMacro.h"
#include "Framework/FrameworkMacro.h"

// System includes
//

// External includes
//

// Class include
//
#include "Simulation/Simulation.hpp"

// Project includes
//
#include "Exceptions/Exception.hpp"
#include "IoTools/Formatter.hpp"

namespace GeoMHDiSCC {

   Simulation::Simulation()
      : mExecutionTimer(), mSimRunCtrl()
   {
   }

   Simulation::~Simulation()
   {
   }

   void Simulation::initBase()
   {
      /// \mhdBug Fake implementation

      // Start timer
      this->mExecutionTimer.start();

      // Make sure to catch raised exception in initialisation steps
      try{
         // Initialise the IO system
         this->mSimIoCtrl.init();

         // Initialise the equation parameters
         this->mspEqParams->init(this->mSimIoCtrl.configPhysical());

         // Get number CPU from configuration file
         int nCpu = this->mSimIoCtrl.configNCpu();

         // Initialise the workflow
         FrameworkMacro::setup(nCpu);
      }
      catch(Exception &e)
      {
         std::cout << e.what() << std::endl;

         throw -1;
      }

      // Stop timer and update initialisation time
      this->mExecutionTimer.stop();
      this->mExecutionTimer.update(ExecutionTimer::INIT);

      // Print message to signal successful completion of initialisation step
      if(FrameworkMacro::allowsIO())
      {
         IoTools::Formatter::printLine(std::cout, '-');
         IoTools::Formatter::printCentered(std::cout, "Simulation initialisation successfull", '*');
         IoTools::Formatter::printLine(std::cout, '-');
         IoTools::Formatter::printNewline(std::cout);
      }

      // Make sure nodes are synchronised after initialisation
      FrameworkMacro::synchronize();
   }

   void Simulation::init()
   {
      // Initialise the variables
      this->initVariables();

      // Setup the equations
      this->setupEquations();

      // Initialise the timestepper
      this->initTimestepper();

      // Setup output files (ASCII diagnostics, state files, etc)
      this->setupOutput();

      // Cleanup IO control
      this->mSimIoCtrl.cleanup();
   }

   void Simulation::run()
   {
      // Start timer
      this->mExecutionTimer.start();

      // Execute pre-run steps
      this->preRun();

      // Stop pre-run timing
      this->mExecutionTimer.stop();
      this->mExecutionTimer.update(ExecutionTimer::PRERUN);

      // Start timer
      this->mExecutionTimer.start();

      // Start main loop of simulation
      while(this->mSimRunCtrl.status() == Runtime::Status::GOON)
      {
         // Compute the nonlinear terms
         this->computeNonlinear();

         // Timestep the equations
         this->timestepEquations();

         // Write the output
         this->writeOutput();

         // Synchronise computation nodes
         FrameworkMacro::synchronize();
      }

      // Stop main loop timing
      this->mExecutionTimer.stop();
      this->mExecutionTimer.update(ExecutionTimer::RUN);

      // Start timer
      this->mExecutionTimer.start();

      // Execute post-run operations
      this->postRun();

      // Synchronise computation nodes
      FrameworkMacro::synchronize();

      // Stop post-run timing
      this->mExecutionTimer.stop();
      this->mExecutionTimer.update(ExecutionTimer::POSTRUN);

      // Synchronise computation nodes
      FrameworkMacro::synchronize();
   }

   void Simulation::finalize()
   {
      // Print execution timer infos
      this->mExecutionTimer.printInfo(std::cout);

      // Print profiling infos (if required)
      ProfilerMacro_printInfo();

      // Print storage profiling infos (if required)
      StorageProfilerMacro_printInfo();
   }

   void Simulation::setInitialStateFile(int spInitFile)//SharedStateFile spInitFile)
   {
      /// \mhdBug Fake implementation
   }

   void Simulation::addOutputFile(int spOutFile)//SharedAscii spOutFile)
   {
      /// \mhdBug Fake implementation
   }

   void Simulation::addOutputFile(double spOutFile)//SharedHdf5 spOutFile)
   {
      /// \mhdBug Fake implementation
   }

   void Simulation::preRun()
   {
      // Print message to signal successful completion of initialisation step
      if(FrameworkMacro::allowsIO())
      {
         IoTools::Formatter::printLine(std::cout, '-');
         IoTools::Formatter::printCentered(std::cout, "... Starting simulation ...", '*');
         IoTools::Formatter::printLine(std::cout, '-');
         IoTools::Formatter::printNewline(std::cout);
      }

      // Write initial ASCII output
      this->mSimIoCtrl.writeAscii();

      // Write initial state file
      this->mSimIoCtrl.writeHdf5();

      // Synchronise all nodes of simulation
      FrameworkMacro::synchronize();
   }

   void Simulation::computeNonlinear()
   {
      /// \mhdBug Fake implementation
      
      // Compute backward transform
      ProfilerMacro_start(ProfilerMacro::BWDTRANSFORM);
      this->mspBwdGrouper->transform(this->mScalarVariables, this->mVectorVariables, this->mTransformCoordinator);
      ProfilerMacro_stop(ProfilerMacro::BWDTRANSFORM);

      // compute nonlinear interaction and forward transform
      this->mspFwdGrouper->transform(this->mScalarEquations, this->mVectorEquations, this->mTransformCoordinator);
   }

   void Simulation::timestepEquations()
   {
      /// \mhdBug Fake implementation
      
      ProfilerMacro_start(ProfilerMacro::TIMESTEP);
      //this->mTimestepper.stepForward(this->mScalarEquations, this->mVectorEquations);
      ProfilerMacro_stop(ProfilerMacro::TIMESTEP);

      ProfilerMacro_start(ProfilerMacro::CONTROL);
      if(this->mTimestepper.finishedStep())
      {
         // Update timestepper
         this->mTimestepper.update();
      
         // Update simulation control
         this->mSimRunCtrl.update();
      }
      ProfilerMacro_stop(ProfilerMacro::CONTROL);
   }

   void Simulation::writeOutput()
   {
      ProfilerMacro_start(ProfilerMacro::IO);
      if(this->mTimestepper.finishedStep() && this->mSimRunCtrl.doIO())
      {
         // Write initial ASCII output
         this->mSimIoCtrl.writeAscii();
      
         // Write initial state file
         this->mSimIoCtrl.writeHdf5();
      }
      ProfilerMacro_stop(ProfilerMacro::IO);
   }

   void Simulation::postRun()
   {
      // Write final ASCII output
      this->mSimIoCtrl.writeAscii();

      // Write final state file
      this->mSimIoCtrl.writeHdf5();

      // Synchronise all nodes of simulation
      FrameworkMacro::synchronize();
   }

   void Simulation::initVariables()
   {
      // Iterator over info
      VariableRequirement::const_iterator infoIt;

      // Storage for the variable info
      VariableRequirement varInfo;

      //
      // Identify the required variables
      //

      // Loop over all scalar equations
      std::vector<SharedIScalarEquation>::iterator scalEqIt;
      for(scalEqIt = this->mScalarEquations.begin(); scalEqIt < this->mScalarEquations.end(); scalEqIt++)
      {
         varInfo.merge((*scalEqIt)->requirements());
      }

      // Loop over all vector equations
      std::vector<SharedIVectorEquation>::iterator vectEqIt;
      for(vectEqIt = this->mVectorEquations.begin(); vectEqIt < this->mVectorEquations.end(); vectEqIt++)
      {
         varInfo.merge((*vectEqIt)->requirements());
      }

      // 
      // Create the required variables
      //

      // Initialise variables
      for(infoIt = varInfo.begin(); infoIt != varInfo.end(); infoIt++)
      {
         // Check if spectral variable is required
         if(infoIt->second.needSpectral())
         {
            // Separate scalar and vector fields
            if(infoIt->second.isScalar())
            {
               // Create the shared scalar variable
               this->mScalarVariables.insert(std::make_pair(infoIt->first, Datatypes::SharedScalarVariableType(new Datatypes::ScalarVariableType(this->mspRes))));
            } else
            {
               // Create the shared vector variable
               this->mVectorVariables.insert(std::make_pair(infoIt->first, Datatypes::SharedVectorVariableType(new Datatypes::VectorVariableType(this->mspRes))));
            }

            // Initialise the physical values if required
            if(infoIt->second.needPhysical())
            {
               // Separate scalar and vector fields
               if(infoIt->second.isScalar())
               {
                  this->mScalarVariables.at(infoIt->first)->initializePhysical();
               } else
               {
                  this->mVectorVariables.at(infoIt->first)->initializePhysical();
               }
            }

            // Initialise the physical differential values if required (gradient or curl)
            if(infoIt->second.needPhysicalDiff())
            {
               // Separate scalar and vector fields
               if(infoIt->second.isScalar())
               {
                  this->mScalarVariables.at(infoIt->first)->initializePhysicalDiff();
               } else
               {
                  this->mVectorVariables.at(infoIt->first)->initializePhysicalDiff();
               }
            }

            // Separate scalar and vector fields
            if(infoIt->second.isScalar())
            {
               // Initialise to zero
               this->mScalarVariables.at(infoIt->first)->setZeros();

               #ifdef EPMPHOENIX_STORAGEPROFILE
                  StorageProfilerMacro_update(StorageProfiler::VARIABLES, this->mScalarVariables.at(infoIt->first)->requiredStorage());
               #endif // EPMPHOENIX_STORAGEPROFILE
            } else
            {
               // Initialise to zero
               this->mVectorVariables.at(infoIt->first)->setZeros();

               #ifdef EPMPHOENIX_STORAGEPROFILE
                  StorageProfilerMacro_update(StorageProfiler::VARIABLES, this->mVectorVariables.at(infoIt->first)->requiredStorage());
               #endif // EPMPHOENIX_STORAGEPROFILE
            }
         }
      }
      
      // Initialise the transform coordinator
      this->mTransformCoordinator.initTransforms(this->mspRes, varInfo);

      // Initialise the communicator
      this->mTransformCoordinator.initCommunicator(this->mspRes);

      // Get the buffer pack sizes
      ArrayI packs1DFwd = this->mspFwdGrouper->packs1D(varInfo);
      ArrayI packs2DFwd = this->mspFwdGrouper->packs2D(varInfo);
      ArrayI packs1DBwd = this->mspBwdGrouper->packs1D(varInfo);
      ArrayI packs2DBwd = this->mspBwdGrouper->packs2D(varInfo);

      // Initialise the converters
      this->mTransformCoordinator.communicator().initConverter(this->mspRes, packs1DFwd, packs1DBwd, packs2DFwd, packs2DBwd, this->mspFwdGrouper->split);
   }

   void Simulation::setupEquations()
   {
      // Loop over all scalar variables
      std::map<PhysicalNames::Id, Datatypes::SharedScalarVariableType>::iterator scalIt;
      for(scalIt = this->mScalarVariables.begin(); scalIt != this->mScalarVariables.end(); scalIt++)
      {
         // Loop over scalar equations
         std::vector<SharedIScalarEquation>::iterator scalEqIt;
         for(scalEqIt = this->mScalarEquations.begin(); scalEqIt < this->mScalarEquations.end(); scalEqIt++)
         {
            // Set scalar variable as unknown scalar field
            if((*scalEqIt)->name() == scalIt->first)
            {
               (*scalEqIt)->setUnknown(this->mScalarVariables.at(scalIt->first));

               // Finish initialisation of equation
               (*scalEqIt)->init();
            }

            // Set scalar variable as additional scalar field
            if((*scalEqIt)->requirements(scalIt->first).needPhysical() || (*scalEqIt)->requirements(scalIt->first).needPhysicalDiff())
            {
               (*scalEqIt)->setField(scalIt->first, this->mScalarVariables.at(scalIt->first));
            }
         }

         // Loop over vector equations
         std::vector<SharedIVectorEquation>::iterator vectEqIt;
         for(vectEqIt = this->mVectorEquations.begin(); vectEqIt < this->mVectorEquations.end(); vectEqIt++)
         {
            // Set scalar variable as additional scalar field
            if((*vectEqIt)->requirements(scalIt->first).needPhysical() || (*vectEqIt)->requirements(scalIt->first).needPhysicalDiff())
            {
               (*vectEqIt)->setField(scalIt->first, this->mScalarVariables.at(scalIt->first));
            }
         }
      }

      // Loop over all vector variables
      std::map<PhysicalNames::Id, Datatypes::SharedVectorVariableType>::iterator vectIt;
      for(vectIt = this->mVectorVariables.begin(); vectIt != this->mVectorVariables.end(); vectIt++)
      {
         // Loop over scalar equations
         std::vector<SharedIScalarEquation>::iterator scalEqIt;
         for(scalEqIt = this->mScalarEquations.begin(); scalEqIt < this->mScalarEquations.end(); scalEqIt++)
         {
            // Set vector variable as additional vector field
            if((*scalEqIt)->requirements(vectIt->first).needPhysical() || (*scalEqIt)->requirements(vectIt->first).needPhysicalDiff())
            {
               (*scalEqIt)->setField(vectIt->first, this->mVectorVariables.at(vectIt->first));
            }
         }

         // Loop over vector equations
         std::vector<SharedIVectorEquation>::iterator vectEqIt;
         for(vectEqIt = this->mVectorEquations.begin(); vectEqIt < this->mVectorEquations.end(); vectEqIt++)
         {
            // Set vector variable as unknown vector field
            if((*vectEqIt)->name() == vectIt->first)
            {
               (*vectEqIt)->setUnknown(this->mVectorVariables.at(vectIt->first));

               // Finish initialisation of equation
               (*vectEqIt)->init();
            }

            // Set vector variable as additional vector field
            if((*vectEqIt)->requirements(vectIt->first).needPhysical() || (*vectEqIt)->requirements(vectIt->first).needPhysicalDiff())
            {
               (*vectEqIt)->setField(vectIt->first, this->mVectorVariables.at(vectIt->first));
            }
         }
      }
   }

   void Simulation::initTimestepper()
   {
//      std::vector<SharedIScalarEquation>::iterator scalEqIt;
//      std::vector<SharedIVectorEquation>::iterator vectEqIt;
//
//      // Storage for the dimension and parameters
//      ArrayI dims(3);
//      dims.setConstant(1);
//      ArrayI specIdx(3);
//      specIdx.setConstant(-1);
//
//      // Create spectral operators for the three dimensions
//      Code::SpectralOperator1DType  spec1D(2, dims(0), specIdx);
//      Code::SpectralOperator2DType  spec2D(2, dims(2), specIdx);
//      Code::SpectralOperator3DType  spec3D(2, dims(1), specIdx);
//
//      bool status = false;
//      bool newMatrix = false;
//
//      dims(2) = this->mspRes->cpu()->dim(0)->dim3D();
//      for(int k = 0; k < dims(2); k++)
//      {
//         specIdx(2) = this->mspRes->cpu()->dim(0)->idx3D(k); 
//         dims(1) = this->mspRes->cpu()->dim(0)->dim2D(k);
//         for(int j = 0; j < dims(1); j++)
//         {
//            specIdx(1) = this->mspRes->cpu()->dim(0)->idx2D(j,k); 
//            dims(0) = this->mspRes->cpu()->dim(0)->dimBwd(j,k);
//            for(int i = 0; i < dims(0); i++)
//            {
//               specIdx(0) = this->mspRes->cpu()->dim(0)->idxBwd(i,j,k); 
//
//               // Reset spectral operator 1D
//               status = spec1D.loopNext(2*dims(0)/3, specIdx);
//               newMatrix = status;
//               // Reset spectral operator 2D
//               status = spec2D.loopNext(dims(2), specIdx);
//               newMatrix = newMatrix || status;
//               // Reset spectral operator 3D
//               status = spec3D.loopNext(dims(1), specIdx);
//               newMatrix = newMatrix || status;
//               if(newMatrix)
//               {
//                  // Loop over all scalar equations
//                  for(scalEqIt = this->mScalarEquations.begin(); scalEqIt < this->mScalarEquations.end(); scalEqIt++)
//                  {
//                     (*scalEqIt)->setSpectralMatrices(spec1D, spec2D, spec3D);
//                  }
//
//                  // Loop over all vector equations
//                  for(vectEqIt = this->mVectorEquations.begin(); vectEqIt < this->mVectorEquations.end(); vectEqIt++)
//                  {
//                     (*vectEqIt)->setSpectralMatrices(spec1D, spec2D, spec3D);
//                  }
//               }
//            }
//         }
//      }
//
//      // Loop over all scalar equations
//      for(scalEqIt = this->mScalarEquations.begin(); scalEqIt < this->mScalarEquations.end(); scalEqIt++)
//      {
//         (*scalEqIt)->finalizeMatrices();
//      }
//
//      // Loop over all vector equations
//      for(vectEqIt = this->mVectorEquations.begin(); vectEqIt < this->mVectorEquations.end(); vectEqIt++)
//      {
//         (*vectEqIt)->finalizeMatrices();
//      }
//
//      // Init timestepper
//      this->timestepper().init(this->mScalarEquations, this->mVectorEquations);
   }

   void Simulation::setupOutput()
   {
//      // Loop over all scalar equations to add scalars
//      std::vector<SharedIScalarEquation>::iterator scalEqIt;
//      for(scalEqIt = this->mScalarEquations.begin(); scalEqIt < this->mScalarEquations.end(); scalEqIt++)
//      {
//         pOutState->addScalar((*this->mScalarVariables.find((*scalEqIt)->name())));
//         if((*scalEqIt)->requirements((*scalEqIt)->name())(1))
//         {
//            pVisState->addScalar((*this->mScalarVariables.find((*scalEqIt)->name())));
//         }
//      }
//
//      // Loop over all vector equations to add vectors
//      std::vector<SharedIVectorEquation>::iterator vectEqIt;
//      for(vectEqIt = this->mVectorEquations.begin(); vectEqIt < this->mVectorEquations.end(); vectEqIt++)
//      {
//         pOutState->addVector((*this->mVectorVariables.find((*vectEqIt)->name())));
//         if((*vectEqIt)->requirements((*vectEqIt)->name())(1))
//         {
//            pVisState->addVector((*this->mVectorVariables.find((*vectEqIt)->name())));
//         }
//      }
//
//      // init the output writers
//      this->mSimIoCtrl.initWriters();
   }

}
