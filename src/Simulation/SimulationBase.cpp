/** \file SimulationBase.cpp
 *  \brief Source of the building block for the implementation of a simulation
 */

// Debug includes
//
#include "Debug/PrepMacros/ProfilerMacro.h"
#include "Debug/PrepMacros/StorageProfilerMacro.h"

// Configuration includes
//
#include "Base/PrepMacros/SmartPointerMacro.h"
#include "Base/PrepMacros/FrameworkMacro.h"
#include "Simulation/PrepMacros/SpatialSchemeMacro.h"
#include "Simulation/PrepMacros/EquationParametersMacro.h"
#include "Simulation/PrepMacros/SpectralOperatorTypedefsMacro.h"

// System includes
//

// External includes
//

// Class include
//
#include "Simulation/System/SimulationBase.hpp"

// Project includes
//
#include "Base/Resolutions/Resolution.hpp"
#include "Base/IO/Ascii/FormatToolbox.hpp"
#include "Base/LoadSplitter/LoadSplitter.hpp"
#include "Base/IO/Control/ConfigParts/PhysicalPart.hpp"
#include "Simulation/IO/Hdf5/StateFileReader.hpp"
#include "Base/LoadSplitter/Algorithms/SplittingDescription.hpp"


namespace EPMPhoenix {

   SimulationBase::SimulationBase()
      : SpatialSimulationBase()
   {
   }

   void SimulationBase::initSimulation()
   {
      // Create the equation parameter shared pointer
      SharedEquationParameters spEqParams(new EquationParametersMacro());
      this->mspEqParams = spEqParams;

      // Add the equation parameter dependent configuration to file
      SharedPhysicalPart   spPhys(new PhysicalPart(this->spEqParams()->names()));
      this->ioSys().spCfg()->addPart(SimulationBlocks::PHYSICAL, spPhys);

      // Initialise the IO system
      this->ioSys().init();

      // Initialise the equation parameters
      this->spEqParams()->init(this->ioSys().spCfg()->physical()->fMap());

      // Get number CPU from configuration file
      int nCpu = this->ioSys().spCfg()->parallel()->iValue("cpus");

      // Extract dimensions from configuration file
      std::map<std::string,int>  trunc = this->ioSys().spCfg()->truncation()->iMap();
      ArrayI dim(trunc.size());
      std::map<std::string,int>::const_iterator  itI;
      int i = 0;
      for(itI = trunc.begin(); itI != trunc.end(); itI++)
      {
         dim(i) = itI->second;
         i++;
      }

      // Initialise spatial component
      this->initSpatial(nCpu, dim);

      // Initialise the simulation control
      this->simCtrl().init();
   }

   void SimulationBase::initTimestepper()
   {
      std::vector<SharedScalarEquation>::iterator scalEqIt;
      std::vector<SharedVectorEquation>::iterator vectEqIt;

      // Storage for the dimension and parameters
      ArrayI dims(3);
      dims.setConstant(1);
      ArrayI specIdx(3);
      specIdx.setConstant(-1);

      // Create spectral operators for the three dimensions
      Code::SpectralOperator1DType  spec1D(2, dims(0), specIdx);
      Code::SpectralOperator2DType  spec2D(2, dims(2), specIdx);
      Code::SpectralOperator3DType  spec3D(2, dims(1), specIdx);

      bool status = false;
      bool newMatrix = false;

      dims(2) = this->spRes()->cpu()->dim(0)->dim3D();
      for(int k = 0; k < dims(2); k++)
      {
         specIdx(2) = this->spRes()->cpu()->dim(0)->idx3D(k); 
         dims(1) = this->spRes()->cpu()->dim(0)->dim2D(k);
         for(int j = 0; j < dims(1); j++)
         {
            specIdx(1) = this->spRes()->cpu()->dim(0)->idx2D(j,k); 
            dims(0) = this->spRes()->cpu()->dim(0)->dimBwd(j,k);
            for(int i = 0; i < dims(0); i++)
            {
               specIdx(0) = this->spRes()->cpu()->dim(0)->idxBwd(i,j,k); 

               // Reset spectral operator 1D
               status = spec1D.loopNext(2*dims(0)/3, specIdx);
               newMatrix = status;
               // Reset spectral operator 2D
               status = spec2D.loopNext(dims(2), specIdx);
               newMatrix = newMatrix || status;
               // Reset spectral operator 3D
               status = spec3D.loopNext(dims(1), specIdx);
               newMatrix = newMatrix || status;
               if(newMatrix)
               {
                  // Loop over all scalar equations
                  for(scalEqIt = this->mScalarEquations.begin(); scalEqIt < this->mScalarEquations.end(); scalEqIt++)
                  {
                     (*scalEqIt)->setSpectralMatrices(spec1D, spec2D, spec3D);
                  }

                  // Loop over all vector equations
                  for(vectEqIt = this->mVectorEquations.begin(); vectEqIt < this->mVectorEquations.end(); vectEqIt++)
                  {
                     (*vectEqIt)->setSpectralMatrices(spec1D, spec2D, spec3D);
                  }
               }
            }
         }
      }

      // Loop over all scalar equations
      for(scalEqIt = this->mScalarEquations.begin(); scalEqIt < this->mScalarEquations.end(); scalEqIt++)
      {
         (*scalEqIt)->finaliseMatrices();
      }

      // Loop over all vector equations
      for(vectEqIt = this->mVectorEquations.begin(); vectEqIt < this->mVectorEquations.end(); vectEqIt++)
      {
         (*vectEqIt)->finaliseMatrices();
      }

      // Init timestepper
      this->timestepper().init(this->mScalarEquations, this->mVectorEquations);
   }

   void SimulationBase::initInput()
   {
//      // Create a state file reader
//      StateFileReader  inState(this->spRes(), "_initial");
//
//      // Add scalars
//      std::map<PhysicalNames::Id,Code::SharedScalarVariable>::const_iterator  sit;
//      for(sit = this->mScalarVariables.begin(); sit != this->mScalarVariables.end(); sit++)
//      {
//         inState.addScalar(*sit);
//      }
//
//      // Add vectors
//      std::map<PhysicalNames::Id,Code::SharedVectorVariable>::const_iterator  vit;
//      for(vit = this->mVectorVariables.begin(); vit != this->mVectorVariables.end(); vit++)
//      {
//         inState.addVector(*vit);
//      }
//
//      // Initialise file
//      inState.init();
//
//      // Read in data
//      inState.read();
//
//      // Finalise file
//      inState.finalise();

      std::map<PhysicalNames::Id,Code::SharedScalarVariable>::const_iterator  sit;
      for(sit = this->mScalarVariables.begin(); sit != this->mScalarVariables.end(); sit++)
      {
         for(int i = 1; i < sit->second->dom(0).spRes()->cpu()->dim(0)->dim3D(); i++)
         {
            sit->second->rDom(0).rPerturbation().rSlice(i).setZero();
            sit->second->rDom(0).rPerturbation().rSlice(i).block(0,0,5,5).setRandom();
            sit->second->rDom(0).rPerturbation().rSlice(i) *= 1.0e-12;
         }
         //sit->second->rDom(0).rPerturbation().rSlice(1)(1,1) = 1.0e-4;
         //sit->second->rDom(0).rPerturbation().rSlice(1)(2,0) = -2*1.0e-4;
         //sit->second->rDom(0).rPerturbation().rSlice(1)(3,1) = -1.0e-4;
         //sit->second->rDom(0).rPerturbation().rSlice(1)(4,0) = 2*1.0e-4;
      }
   }

   void SimulationBase::initBaseOutput()
   {
   }

   void SimulationBase::setupOutput()
   {
      // init the output writers
      this->mIOSystem.initWriters();
   }

   void SimulationBase::cleanupSimulation()
   {
      // Initialise the base system
      this->cleanupSystem();

      // Initialise the IO system
      this->mIOSystem.cleanup();
   }

   void SimulationBase::finaliseOutput()
   {
      // finalise the output writers
      this->mIOSystem.finaliseWriters();
   }

   void SimulationBase::finaliseSimulation()
   {
      // Finalise IOSystem
      this->mIOSystem.finalise();
   }

   void SimulationBase::timestepEquations()
   {
      ProfilerMacro_start(ProfilerMacro::TIMESTEP);
      this->timestepper().stepForward(this->mScalarEquations, this->mVectorEquations);
      ProfilerMacro_stop(ProfilerMacro::TIMESTEP);

      ProfilerMacro_start(ProfilerMacro::CONTROL);
      if(this->timestepper().finishedStep())
      {
         // Update timestepper
         this->timestepper().update();
      
         // Update simulation control
         this->simCtrl().update();
      }
      ProfilerMacro_stop(ProfilerMacro::CONTROL);
   }

   void SimulationBase::writeOutput()
   {
      ProfilerMacro_start(ProfilerMacro::IO);
      if(this->timestepper().finishedStep() && this->simCtrl().doIO())
      {
         // Write initial ASCII output
         this->ioSys().writeASCII();

         // Write initial state file
         this->ioSys().writeHdf5();
      }
      ProfilerMacro_stop(ProfilerMacro::IO);
   }

   void SimulationBase::preRun()
   {
      // Print message to signal successful completion of initialisation step
      if(FrameworkMacro::allowsIO())
      {
         FormatToolbox::printLine('-');
         FormatToolbox::printCentered("... Starting simulation ...", '*');
         FormatToolbox::printLine('-');
         FormatToolbox::printNewline();
      }

      // Write initial ASCII output
      this->ioSys().writeASCII();

      // Write initial state file
      this->ioSys().writeHdf5();

      // Synchronise all nodes of simulation
      FrameworkMacro::synchronise();

      // Execute last initialisations
      this->simCtrl().preRun();
   }

   void SimulationBase::postRun()
   {
      // Write final ASCII output
      this->ioSys().writeASCII();

      // Write final state file
      this->ioSys().writeHdf5();

      // Synchronise all nodes of simulation
      FrameworkMacro::synchronise();

      // Execute run finalisations
      this->simCtrl().postRun();
   }

}
