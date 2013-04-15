/** \file Timestepper.cpp
 *  \brief Implementation of a general timestepper structure
 */

// Debug includes
//

// System includes
//

// External includes
//
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>

// Class include
//
#include "Timesteppers/Timestepper.hpp"

// Project includes
//
#include "Base/MathConstants.hpp"
#include "Timesteppers/ImExRK3.hpp"

#include "Timers/TimerMacro.h"
#include <iostream>
namespace GeoMHDiSCC {

namespace Timestep {

   Timestepper::Timestepper()
      : mcMaxJump(1.602), mcUpWindow(1.05), mcMinDt(1e-8), mStep(0), mOldDt(this->mcMinDt), mDt(this->mcMinDt), mTime(0.0)
   {
   }

   Timestepper::~Timestepper()
   {
   }

   MHDFloat Timestepper::time() const
   {
      return this->mTime;
   }

   MHDFloat Timestepper::timestep() const
   {
      return this->mDt;
   }

   bool Timestepper::finishedStep() const
   {
      return this->mStep == 0;
   }

   void Timestepper::update()
   {
      this->mTime += this->mDt;
   }

   void Timestepper::adaptTimestep(const MHDFloat cfl, const std::vector<Equations::SharedIScalarEquation>& scalEq, const std::vector<Equations::SharedIVectorEquation>& vectEq)
   {
      // Flag to update timestep
      bool hasNewDt = false;

      // Check if CFL allows for a larger timestep
      if(cfl > this->mcUpWindow*this->mDt)
      {
         // Activate matrices update
         hasNewDt = true;

         // Set new timestep
         this->mOldDt = this->mDt;
         this->mDt = std::min(cfl, this->mcMaxJump*this->mDt);
      
      // Check if CFL is below minimal timestep or downard jump is large
      } else if(cfl < this->mcMinDt || cfl < this->mDt/this->mcMaxJump)
      {
         // Don't update matrices
         hasNewDt = false;
 
         // Signal simulation abort
         this->mOldDt = this->mDt;
         this->mDt = -cfl;
     
      // Check if CFL requires a lower timestep
      } else if(cfl < this->mDt)
      {
         // Activate matrices update
         hasNewDt = true;

         // Set new timestep
         this->mOldDt = this->mDt;
         this->mDt = cfl/this->mcUpWindow;

      // No need to change timestep
      } else
      {
         hasNewDt = false;
      }

      if(hasNewDt)
      {
         std::cerr << "Updating timestep and matrices with new Dt = " << this->mDt << std::endl;
      }

      //
      // Update the timestep matrices if necessary
      //
      if(hasNewDt)
      {
TimerMacro  timer;
timer.start();
std::cerr << "START update matrices" << std::endl;
         // Loop over all substeps of timestepper
         for(this->mStep = 0; this->mStep < ImExRK3::STEPS; this->mStep++)
         {
            // Loop over all scalar equations
            std::vector<Equations::SharedIScalarEquation>::const_iterator scalEqIt;
            for(scalEqIt = scalEq.begin(); scalEqIt < scalEq.end(); scalEqIt++)
            {
               // Create (coupled) matrices
               this->updateMatrices((*scalEqIt), FieldComponents::Spectral::SCALAR);
            }

            // Loop over all vector equations
            std::vector<Equations::SharedIVectorEquation>::const_iterator vectEqIt;
            for(vectEqIt = vectEq.begin(); vectEqIt < vectEq.end(); vectEqIt++)
            {
               // Create (coupled) matrices
               this->updateMatrices((*vectEqIt), FieldComponents::Spectral::ONE);

               // Create (coupled) matrices
               this->updateMatrices((*vectEqIt), FieldComponents::Spectral::TWO);
            }
         }
timer.stop();
std::cerr << "FINISH update matrices: " << timer.time() << std::endl;

timer.start();
std::cerr << "START complex solver update" << std::endl;
         // Update solvers from complex equation steppers
         for(size_t i = 0; i < this->mEqZStepper.size(); i++)
         {
            this->mEqZStepper.at(i).updateSolver();
         }
timer.stop();
std::cerr << "FINISH complex solver update: " << timer.time() << std::endl;

timer.start();
std::cerr << "START real solver update" << std::endl;
         // Update solvers from real equation steppers
         for(size_t i = 0; i < this->mEqDStepper.size(); i++)
         {
            this->mEqDStepper.at(i).updateSolver();
         }
timer.stop();
std::cerr << "FINISH real solver update: " << timer.time() << std::endl;

         // Reset the step index
         this->mStep = 0;
      }
   }

   void Timestepper::stepForward(const std::vector<Equations::SharedIScalarEquation>& scalEq, const std::vector<Equations::SharedIVectorEquation>& vectEq)
   {
      // Update the equation input to the timestepper
      this->getInput(scalEq, vectEq);

      // Compute the RHS of the linear systems
      this->computeRHS();

      // Solve all the linear systems
      this->solve();
      
      // Transfer timestep output back to equations
      this->transferOutput(scalEq, vectEq);

      // Update the internal step counter, counting from 0 to steps - 1
      this->mStep = (this->mStep + 1) % ImExRK3::STEPS;
   }

   void Timestepper::init(const MHDFloat dt, const std::vector<Equations::SharedIScalarEquation>& scalEq, const std::vector<Equations::SharedIVectorEquation>& vectEq)
   {
      // Set initial timestep
      this->mOldDt = dt;
      this->mDt = dt;

      //
      // Determine if solvers need to be complex or not
      //

      // Loop over all scalar equations
      std::vector<Equations::SharedIScalarEquation>::const_iterator scalEqIt;
      for(scalEqIt = scalEq.begin(); scalEqIt < scalEq.end(); scalEqIt++)
      {
         // Get type information for the equation steppers
         this->getEqStepperType((*scalEqIt), FieldComponents::Spectral::SCALAR);
      }

      // Loop over all vector equations
      std::vector<Equations::SharedIVectorEquation>::const_iterator vectEqIt;
      for(vectEqIt = vectEq.begin(); vectEqIt < vectEq.end(); vectEqIt++)
      {
         // Get type information for the equation steppers for the toroidal component
         this->getEqStepperType((*vectEqIt), FieldComponents::Spectral::ONE);

         // Get type information for the equation steppers for the poloidal component
         this->getEqStepperType((*vectEqIt), FieldComponents::Spectral::TWO);
      }

      //
      // Create real/complex timesteppers and update coupling information
      //

      // Loop over all scalar equations
      for(scalEqIt = scalEq.begin(); scalEqIt < scalEq.end(); scalEqIt++)
      {
         // Get type information for the equation steppers
         this->createEqStepper((*scalEqIt), FieldComponents::Spectral::SCALAR);
      }

      // Loop over all vector equations
      for(vectEqIt = vectEq.begin(); vectEqIt < vectEq.end(); vectEqIt++)
      {
         // Get type information for the equation steppers for the toroidal component
         this->createEqStepper((*vectEqIt), FieldComponents::Spectral::ONE);

         // Get type information for the equation steppers for the poloidal component
         this->createEqStepper((*vectEqIt), FieldComponents::Spectral::TWO);
      }

      //
      // Create the timestep matrices
      //

TimerMacro  timer;
timer.start();
std::cerr << "START create matrices" << std::endl;
      // Loop over all substeps of timestepper
      for(this->mStep = 0; this->mStep < ImExRK3::STEPS; this->mStep++)
      {
         // Loop over all scalar equations
         for(scalEqIt = scalEq.begin(); scalEqIt < scalEq.end(); scalEqIt++)
         {
            // Create (coupled) matrices
            this->createMatrices((*scalEqIt), FieldComponents::Spectral::SCALAR);
         }

         // Loop over all vector equations
         for(vectEqIt = vectEq.begin(); vectEqIt < vectEq.end(); vectEqIt++)
         {
            // Create (coupled) matrices
            this->createMatrices((*vectEqIt), FieldComponents::Spectral::ONE);

            // Create (coupled) matrices
            this->createMatrices((*vectEqIt), FieldComponents::Spectral::TWO);
         }
      }
timer.stop();
std::cerr << "FINISH create matrices: " << timer.time() << std::endl;

      //
      // Initialise the solvers and the initial state
      //

timer.start();
std::cerr << "START complex solver init" << std::endl;
      // Initialise solvers from complex equation steppers
      for(size_t i = 0; i < this->mEqZStepper.size(); i++)
      {
         this->mEqZStepper.at(i).initSolver();
      }
timer.stop();
std::cerr << "FINISH complex solver init: " << timer.time() << std::endl;

timer.start();
std::cerr << "START real solver init" << std::endl;
      // Initialise solvers from real equation steppers
      for(size_t i = 0; i < this->mEqDStepper.size(); i++)
      {
         this->mEqDStepper.at(i).initSolver();
      }
timer.stop();
std::cerr << "FINISH real solver init: " << timer.time() << std::endl;

      // Reset the step index
      this->mStep = 0;

      // Initialise with initial state
      this->initSolution(scalEq, vectEq);
   }

   void Timestepper::getEqStepperType(Equations::SharedIEvolutionEquation spEq, FieldComponents::Spectral::Id comp)
   {
      // Get Id for current unknown field
      TimestepCoupling::FieldIdType myId = std::make_pair(spEq->name(),comp);

      // Equation is coupled to other equation and thus already stored
      if(this->mTimeCoupling.isPresent(myId))
      {
         // Get system index for field
         int idx = this->mTimeCoupling.idx(myId);

         // Update complex flag if required
         this->mTimeCoupling.updateType(idx, spEq->isComplex(comp));

         // Check the starting index
         this->mTimeCoupling.checkStart(idx, spEq->startIndex(comp));

      // Equation is not coupled or is the first of coupled set
      } else
      {
         // Get a new index for the system
         int idx = this->mTimeCoupling.newIndex();

         // Add equation unknown to time coupling information
         this->mTimeCoupling.addField(myId, spEq->isComplex(comp), idx, spEq->startIndex(comp));

         // Get the range of coupled fields
         Equations::CouplingInformation::field_iterator it;
         Equations::CouplingInformation::field_iterator_range range = spEq->couplingInfo().fieldRange(comp);

         // Loop over coupled fields and add them to timestep coupling information
         for(it = range.first; it != range.second; ++it)
         {
            this->mTimeCoupling.addField(it->second, false, idx, spEq->startIndex(comp));
         }
      }
   }

   void Timestepper::createEqStepper(Equations::SharedIEvolutionEquation spEq, FieldComponents::Spectral::Id comp)
   {
      // Get Id for current unknown field
      TimestepCoupling::FieldIdType myId = std::make_pair(spEq->name(),comp);

      if(this->mTimeCoupling.idx(myId) < 0)
      {
         // Get the number of coupled fields
         int nC = spEq->couplingInfo().nFields(comp);

         int idx;
         // Create correct equation timestepper
         if(this->mTimeCoupling.isComplex(myId))
         {
            idx = this->mEqZStepper.size();
            this->mEqZStepper.push_back(EquationZTimestepper(nC+1, spEq->startIndex(comp)));
         } else
         {
            idx = this->mEqDStepper.size();
            this->mEqDStepper.push_back(EquationDTimestepper(nC+1, spEq->startIndex(comp)));
         }

         this->mTimeCoupling.updateIndex(this->mTimeCoupling.idx(myId), idx);
      }
   }

   void Timestepper::createMatrices(Equations::SharedIEvolutionEquation spEq, FieldComponents::Spectral::Id comp)
   {
      // ID of the current field
      TimestepCoupling::FieldIdType myId = std::make_pair(spEq->name(),comp);

      // Index of the current field
      int myIdx = this->mTimeCoupling.idx(myId);

      // Get internal information
      std::pair<int,ArrayI>  interInfo = spEq->couplingInfo().internal(comp);

      // Get the number of coupled fields
      int nC = spEq->couplingInfo().nFields(comp);

      // start index for matrices
      int start = this->mStep*interInfo.first;

      // Start row for storage information
      ArrayI startRow(interInfo.first);

      // Storage for the field index
      int fieldIdx;

      // Complex matrices in linear solve
      if(this->mTimeCoupling.isComplex(myId))
      {
         fieldIdx = this->mEqZStepper.at(myIdx).current();

         if(fieldIdx == 0)
         {
            // Reserve space for the matrices to avoid large number of expensive reallocations
            if(this->mStep == 0)
            {
               this->mEqZStepper.at(myIdx).reserveMatrices(ImExRK3::STEPS*interInfo.first);
            }

            int size = 0;
            for(int i = 0; i < interInfo.first; i++)
            {
               // Clear LHS triplets
               this->mEqZStepper.at(myIdx).rLHSTriplets(start+i).clear();
               // Add LHS triplets
               size = this->buildLHSMatrix(this->mEqZStepper.at(myIdx).rLHSTriplets(start+i), spEq, comp, i, nC, fieldIdx);

               // Clear RHS triplets
               this->mEqZStepper.at(myIdx).rRHSTriplets(start+i).clear();
               // Add RHS triplets
               this->buildRHSMatrix(this->mEqZStepper.at(myIdx).rRHSTriplets(start+i), spEq, comp, i, nC, fieldIdx);

               if(this->mStep == 0)
               {
                  // Create RHS and solution data storage
                  this->mEqZStepper.at(myIdx).addStorage(size, interInfo.second(i));
               }

               // Set the start row
               startRow(i) = 0;
            }
         } else
         {
            for(int i = 0; i < interInfo.first; i++)
            {
               // Add LHS triplets
               this->buildLHSMatrix(this->mEqZStepper.at(myIdx).rLHSTriplets(start+i), spEq, comp, i, nC, fieldIdx);

               // Add RHS triplets
               this->buildRHSMatrix(this->mEqZStepper.at(myIdx).rRHSTriplets(start+i), spEq, comp, i, nC, fieldIdx);

               // Set the start row
               startRow(i) = spEq->rowShift(comp, i);
            }
         }

         if(this->mStep == 0)
         {
            // Store storage information
            this->mEqZStepper.at(myIdx).addInformation(myId,startRow);
         }

         // Move field counter to next field
         this->mEqZStepper.at(myIdx).next();

      // Real matrices in linear solve
      } else
      {
         fieldIdx = this->mEqDStepper.at(myIdx).current();

         if(fieldIdx == 0)
         {
            // Reserve space for the matrices to avoid large number of expensive reallocations
            if(this->mStep == 0)
            {
               this->mEqDStepper.at(myIdx).reserveMatrices(ImExRK3::STEPS*interInfo.first);
            }

            int size = 0;
            for(int i = 0; i < interInfo.first; i++)
            {
               // Clear LHS triplets
               this->mEqDStepper.at(myIdx).rLHSTriplets(start+i).clear();
               // Add LHS triplets
               size = this->buildLHSMatrix(this->mEqDStepper.at(myIdx).rLHSTriplets(start+i), spEq, comp, i, nC, fieldIdx);

               // Clear RHS triplets
               this->mEqDStepper.at(myIdx).rRHSTriplets(start+i).clear();
               // Add RHS triplets
               this->buildRHSMatrix(this->mEqDStepper.at(myIdx).rRHSTriplets(start+i), spEq, comp, i, nC, fieldIdx);

               if(this->mStep == 0)
               {
                  // Create RHS and solution data storage
                  this->mEqDStepper.at(myIdx).addStorage(size, interInfo.second(i));
               }

               // Set the start row
               startRow(i) = 0;
            }
         } else
         {
            for(int i = 0; i < interInfo.first; i++)
            {
               // Add LHS triplets
               this->buildLHSMatrix(this->mEqDStepper.at(myIdx).rLHSTriplets(start+i), spEq, comp, i, nC, fieldIdx);

               // Add RHS triplets
               this->buildRHSMatrix(this->mEqDStepper.at(myIdx).rRHSTriplets(start+i), spEq, comp, i, nC, fieldIdx);

               // Set the start row
               startRow(i) = spEq->rowShift(comp, i);
            }
         }

         if(this->mStep == 0)
         {
            // Store storage information
            this->mEqDStepper.at(myIdx).addInformation(myId,startRow);
         }

         // Move field counter to next field
         this->mEqDStepper.at(myIdx).next();
      }
   }

   void Timestepper::updateMatrices(Equations::SharedIEvolutionEquation spEq, FieldComponents::Spectral::Id comp)
   {
      // ID of the current field
      TimestepCoupling::FieldIdType myId = std::make_pair(spEq->name(),comp);

      // Index of the current field
      int myIdx = this->mTimeCoupling.idx(myId);

      // Get internal information
      std::pair<int,ArrayI>  interInfo = spEq->couplingInfo().internal(comp);

      // Get the number of coupled fields
      int nC = spEq->couplingInfo().nFields(comp);

      // start index for matrices
      int start = this->mStep*interInfo.first;

      // Start row for storage information
      ArrayI startRow(interInfo.first);

      // Storage for the field index
      int fieldIdx;

      // Complex matrices in linear solve
      if(this->mTimeCoupling.isComplex(myId))
      {
         fieldIdx = this->mEqZStepper.at(myIdx).current();

         if(fieldIdx == 0)
         {
            for(int i = 0; i < interInfo.first; i++)
            {
               // Update LHS triplets
               this->updateLHSMatrix(this->mEqZStepper.at(myIdx).rLHSTriplets(start+i), spEq, comp, i, nC, fieldIdx);

               // Update RHS triplets
               this->updateRHSMatrix(this->mEqZStepper.at(myIdx).rRHSTriplets(start+i), spEq, comp, i, nC, fieldIdx);

               // Set the start row
               startRow(i) = 0;
            }
         } else
         {
            for(int i = 0; i < interInfo.first; i++)
            {
               // Update LHS triplets
               this->updateLHSMatrix(this->mEqZStepper.at(myIdx).rLHSTriplets(start+i), spEq, comp, i, nC, fieldIdx);

               // Update RHS triplets
               this->updateRHSMatrix(this->mEqZStepper.at(myIdx).rRHSTriplets(start+i), spEq, comp, i, nC, fieldIdx);

               // Set the start row
               startRow(i) = spEq->rowShift(comp, i);
            }
         }

         // Move field counter to next field
         this->mEqZStepper.at(myIdx).next();

      // Real matrices in linear solve
      } else
      {
         fieldIdx = this->mEqDStepper.at(myIdx).current();

         if(fieldIdx == 0)
         {
            for(int i = 0; i < interInfo.first; i++)
            {
               // Update LHS triplets
               this->updateLHSMatrix(this->mEqDStepper.at(myIdx).rLHSTriplets(start+i), spEq, comp, i, nC, fieldIdx);

               // Update RHS triplets
               this->updateRHSMatrix(this->mEqDStepper.at(myIdx).rRHSTriplets(start+i), spEq, comp, i, nC, fieldIdx);

               // Set the start row
               startRow(i) = 0;
            }
         } else
         {
            for(int i = 0; i < interInfo.first; i++)
            {
               // Update LHS triplets
               this->updateLHSMatrix(this->mEqDStepper.at(myIdx).rLHSTriplets(start+i), spEq, comp, i, nC, fieldIdx);

               // Update RHS triplets
               this->updateRHSMatrix(this->mEqDStepper.at(myIdx).rRHSTriplets(start+i), spEq, comp, i, nC, fieldIdx);

               // Set the start row
               startRow(i) = spEq->rowShift(comp, i);
            }
         }

         // Move field counter to next field
         this->mEqDStepper.at(myIdx).next();
      }
   }

   void Timestepper::initSolution(const std::vector<Equations::SharedIScalarEquation>& scalEq, const std::vector<Equations::SharedIVectorEquation>& vectEq)
   {
      // Storage for information and identity
      TimestepCoupling::FieldIdType myId;

      // Loop over all scalar equations
      std::vector<Equations::SharedIScalarEquation>::const_iterator scalEqIt;
      for(scalEqIt = scalEq.begin(); scalEqIt < scalEq.end(); scalEqIt++)
      {
         // Get identity and corresponding equation information
         myId = std::make_pair((*scalEqIt)->name(), FieldComponents::Spectral::SCALAR);

         // Get index of current field
         int myIdx = this->mTimeCoupling.idx(myId);

         // Linear solve matrices are complex
         if(this->mTimeCoupling.isComplex(myId))
         {
            // Get timestep input
            for(int i = 0; i < this->mEqZStepper.at(myIdx).nSystem(); i++)
            {
               (*scalEqIt)->copyTInput(FieldComponents::Spectral::SCALAR, this->mEqZStepper.at(myIdx).rSolution(i), i, this->mEqZStepper.at(myIdx).startRow(myId,i));
            }

         // Linear solve matrices are real
         } else
         {
            // Get timestep input
            for(int i = 0; i < this->mEqDStepper.at(myIdx).nSystem(); i++)
            {
               (*scalEqIt)->copyTInput(FieldComponents::Spectral::SCALAR, this->mEqDStepper.at(myIdx).rSolution(i), i, this->mEqDStepper.at(myIdx).startRow(myId,i));
            }
         }
      }

      // Loop over all vector equations
      std::vector<Equations::SharedIVectorEquation>::const_iterator vectEqIt;
      for(vectEqIt = vectEq.begin(); vectEqIt < vectEq.end(); vectEqIt++)
      {
         // Get identity and corresponding equation information for toroidal component
         myId = std::make_pair((*vectEqIt)->name(), FieldComponents::Spectral::ONE);

         // Get index of current field
         int myIdx = this->mTimeCoupling.idx(myId);

         // Linear solve matrices are complex
         if(this->mTimeCoupling.isComplex(myId))
         {
            // Get timestep input for toroidal component
            for(int i = 0; i < this->mEqZStepper.at(myIdx).nSystem(); i++)
            {
               (*vectEqIt)->copyTInput(FieldComponents::Spectral::ONE, this->mEqZStepper.at(myIdx).rSolution(i), i, this->mEqZStepper.at(myIdx).startRow(myId,i));
            }

         // Linear solve matrices are real
         } else
         {
            // Get timestep input for toroidal component
            for(int i = 0; i < this->mEqDStepper.at(myIdx).nSystem(); i++)
            {
               (*vectEqIt)->copyTInput(FieldComponents::Spectral::ONE, this->mEqDStepper.at(myIdx).rSolution(i), i, this->mEqDStepper.at(myIdx).startRow(myId,i));
            }
         }

         // Get identity and corresponding equation information for second component
         myId = std::make_pair((*vectEqIt)->name(), FieldComponents::Spectral::TWO);

         // Get index of current field
         myIdx = this->mTimeCoupling.idx(myId);

         // Linear solve matrices are complex
         if(this->mTimeCoupling.isComplex(myId))
         {
            // Get timestep input for poloidal component
            for(int i = 0; i < this->mEqZStepper.at(myIdx).nSystem(); i++)
            {
               (*vectEqIt)->copyTInput(FieldComponents::Spectral::TWO, this->mEqZStepper.at(myIdx).rSolution(i), i, this->mEqZStepper.at(myIdx).startRow(myId,i));
            }

         // Linear solve matrices are real
         } else
         {
            // Get timestep input for poloidal component
            for(int i = 0; i < this->mEqDStepper.at(myIdx).nSystem(); i++)
            {
               (*vectEqIt)->copyTInput(FieldComponents::Spectral::TWO, this->mEqDStepper.at(myIdx).rSolution(i), i, this->mEqDStepper.at(myIdx).startRow(myId,i));
            }
         }
      }
   }

   void Timestepper::getInput(const std::vector<Equations::SharedIScalarEquation>& scalEq, const std::vector<Equations::SharedIVectorEquation>& vectEq)
   {
      // Storage for information and identity
      TimestepCoupling::FieldIdType myId;

      // Loop over all scalar equations
      std::vector<Equations::SharedIScalarEquation>::const_iterator scalEqIt;
      for(scalEqIt = scalEq.begin(); scalEqIt < scalEq.end(); scalEqIt++)
      {
         // Get identity and corresponding equation information
         myId = std::make_pair((*scalEqIt)->name(), FieldComponents::Spectral::SCALAR);

         // Get index of current field
         int myIdx = this->mTimeCoupling.idx(myId);

         // Linear solve matrices are complex
         if(this->mTimeCoupling.isComplex(myId))
         {
            // Get timestep input
            for(int i = 0; i < this->mEqZStepper.at(myIdx).nSystem(); i++)
            {
               (*scalEqIt)->timestepInput(FieldComponents::Spectral::SCALAR, this->mEqZStepper.at(myIdx).rRHSData(i), i, this->mEqZStepper.at(myIdx).startRow(myId,i));
            }

         // Linear solve matrices are real
         } else
         {
            // Get timestep input
            for(int i = 0; i < this->mEqDStepper.at(myIdx).nSystem(); i++)
            {
               (*scalEqIt)->timestepInput(FieldComponents::Spectral::SCALAR, this->mEqDStepper.at(myIdx).rRHSData(i), i, this->mEqDStepper.at(myIdx).startRow(myId,i));
            }
         }
      }

      // Loop over all vector equations
      std::vector<Equations::SharedIVectorEquation>::const_iterator vectEqIt;
      for(vectEqIt = vectEq.begin(); vectEqIt < vectEq.end(); vectEqIt++)
      {
         // Get identity and corresponding equation information for toroidal component
         myId = std::make_pair((*vectEqIt)->name(), FieldComponents::Spectral::ONE);

         // Get index of current field
         int myIdx = this->mTimeCoupling.idx(myId);

         // Linear solve matrices are complex
         if(this->mTimeCoupling.isComplex(myId))
         {
            // Get timestep input for toroidal component
            for(int i = 0; i < this->mEqZStepper.at(myIdx).nSystem(); i++)
            {
               (*vectEqIt)->timestepInput(FieldComponents::Spectral::ONE, this->mEqZStepper.at(myIdx).rRHSData(i), i, this->mEqZStepper.at(myIdx).startRow(myId,i));
            }

         // Linear solve matrices are real
         } else
         {
            // Get timestep input for toroidal component
            for(int i = 0; i < this->mEqDStepper.at(myIdx).nSystem(); i++)
            {
               (*vectEqIt)->timestepInput(FieldComponents::Spectral::ONE, this->mEqDStepper.at(myIdx).rRHSData(i), i, this->mEqDStepper.at(myIdx).startRow(myId,i));
            }
         }

         // Get identity and corresponding equation information for poloidal component
         myId = std::make_pair((*vectEqIt)->name(), FieldComponents::Spectral::TWO);

         // Get index of current field
         myIdx = this->mTimeCoupling.idx(myId);

         // Linear solve matrices are complex
         if(this->mTimeCoupling.isComplex(myId))
         {
            // Get timestep input for poloidal component
            for(int i = 0; i < this->mEqZStepper.at(myIdx).nSystem(); i++)
            {
               (*vectEqIt)->timestepInput(FieldComponents::Spectral::TWO, this->mEqZStepper.at(myIdx).rRHSData(i), i, this->mEqZStepper.at(myIdx).startRow(myId,i));
            }

         // Linear solve matrices are real
         } else
         {
            // Get timestep input for poloidal component
            for(int i = 0; i < this->mEqDStepper.at(myIdx).nSystem(); i++)
            {
               (*vectEqIt)->timestepInput(FieldComponents::Spectral::TWO, this->mEqDStepper.at(myIdx).rRHSData(i), i, this->mEqDStepper.at(myIdx).startRow(myId,i));
            }
         }
      }
   }

   void Timestepper::computeRHS()
   {
      // Compute RHS component for complex linear systems
      std::vector<EquationZTimestepper>::iterator   zIt;
      for(zIt = this->mEqZStepper.begin(); zIt != this->mEqZStepper.end(); ++zIt)
      {
         // Compute linear solve RHS
         zIt->computeRHS(this->mStep);
      }

      std::vector<EquationDTimestepper>::iterator   dIt;
      for(dIt = this->mEqDStepper.begin(); dIt != this->mEqDStepper.end(); ++dIt)
      {
         // Compute linear solve RHS
         dIt->computeRHS(this->mStep);
      }
   }

   void Timestepper::solve()
   {
      // Solve complex linear systems
      std::vector<EquationZTimestepper>::iterator   zIt;
      for(zIt = this->mEqZStepper.begin(); zIt != this->mEqZStepper.end(); ++zIt)
      {
         // Compute linear solve RHS
         zIt->solve(this->mStep);
      }

      // Solve real linear systems
      std::vector<EquationDTimestepper>::iterator   dIt;
      for(dIt = this->mEqDStepper.begin(); dIt != this->mEqDStepper.end(); ++dIt)
      {
         // Compute linear solve RHS
         dIt->solve(this->mStep);
      }
   }

   void Timestepper::transferOutput(const std::vector<Equations::SharedIScalarEquation>& scalEq, const std::vector<Equations::SharedIVectorEquation>& vectEq)
   {
      // Storage for information and identity
      TimestepCoupling::FieldIdType myId;

      // Loop over all scalar equations
      std::vector<Equations::SharedIScalarEquation>::const_iterator scalEqIt;
      for(scalEqIt = scalEq.begin(); scalEqIt < scalEq.end(); scalEqIt++)
      {
         // Get identity and corresponding equation information
         myId = std::make_pair((*scalEqIt)->name(), FieldComponents::Spectral::SCALAR);

         // Get index of current field
         int myIdx = this->mTimeCoupling.idx(myId);

         // Linear solve matrices are complex
         if(this->mTimeCoupling.isComplex(myId))
         {
            // Get timestep output
            for(int i = 0; i < this->mEqZStepper.at(myIdx).nSystem(); i++)
            {
               (*scalEqIt)->timestepOutput(FieldComponents::Spectral::SCALAR, this->mEqZStepper.at(myIdx).solution(i), i, this->mEqZStepper.at(myIdx).startRow(myId,i));
            }

         // Linear solve matrices are real
         } else
         {
            // Get timestep output
            for(int i = 0; i < this->mEqDStepper.at(myIdx).nSystem(); i++)
            {
               (*scalEqIt)->timestepOutput(FieldComponents::Spectral::SCALAR, this->mEqDStepper.at(myIdx).solution(i), i, this->mEqDStepper.at(myIdx).startRow(myId,i));
            }
         }
      }

      // Loop over all vector equations
      std::vector<Equations::SharedIVectorEquation>::const_iterator vectEqIt;
      for(vectEqIt = vectEq.begin(); vectEqIt < vectEq.end(); vectEqIt++)
      {
         // Get identity and corresponding equation information for toroidal component
         myId = std::make_pair((*vectEqIt)->name(), FieldComponents::Spectral::ONE);

         // Get index of current field
         int myIdx = this->mTimeCoupling.idx(myId);

         // Linear solve matrices are complex
         if(this->mTimeCoupling.isComplex(myId))
         {
            // Get timestep output for toroidal component
            for(int i = 0; i < this->mEqZStepper.at(myIdx).nSystem(); i++)
            {
               (*vectEqIt)->timestepOutput(FieldComponents::Spectral::ONE, this->mEqZStepper.at(myIdx).solution(i), i, this->mEqZStepper.at(myIdx).startRow(myId,i));
            }

         // Linear solve matrices are real
         } else
         {
            // Get timestep output for toroidal component
            for(int i = 0; i < this->mEqDStepper.at(myIdx).nSystem(); i++)
            {
               (*vectEqIt)->timestepOutput(FieldComponents::Spectral::ONE, this->mEqDStepper.at(myIdx).solution(i), i, this->mEqDStepper.at(myIdx).startRow(myId,i));
            }
         }

         // Get identity and corresponding equation information for poloidal component
         myId = std::make_pair((*vectEqIt)->name(), FieldComponents::Spectral::TWO);

         // Get index of current field
         myIdx = this->mTimeCoupling.idx(myId);

         // Linear solve matrices are complex
         if(this->mTimeCoupling.isComplex(myId))
         {
            // Get timestep input for poloidal component
            for(int i = 0; i < this->mEqZStepper.at(myIdx).nSystem(); i++)
            {
               (*vectEqIt)->timestepOutput(FieldComponents::Spectral::TWO, this->mEqZStepper.at(myIdx).solution(i), i, this->mEqZStepper.at(myIdx).startRow(myId,i));
            }

         // Linear solve matrices are real
         } else
         {
            // Get timestep output for poloidal component
            for(int i = 0; i < this->mEqDStepper.at(myIdx).nSystem(); i++)
            {
               (*vectEqIt)->timestepOutput(FieldComponents::Spectral::TWO, this->mEqDStepper.at(myIdx).solution(i), i, this->mEqDStepper.at(myIdx).startRow(myId,i));
            }
         }
      }
   }

   int Timestepper::buildLHSMatrix(std::vector<Triplet>& triplets, Equations::SharedIEvolutionEquation spEq, FieldComponents::Spectral::Id comp, const int idx, const int nC, const int cRow)
   {
      // Matrix coefficients
      MHDFloat timeCoeff = ImExRK3::lhsT(this->mStep)*1.0/this->mDt;
      MHDFloat linearCoeff = ImExRK3::lhsL(this->mStep);

      // Set block size
      int blockRows = spEq->linearMatrix(comp,idx).first.rows();
      int blockCols = spEq->linearMatrix(comp,idx).first.cols();

      // Block shifts
      int rowShift = 0;
      int colShift = 0;

      // Set block shift for main LHS matrix
      if(nC > 0)
      {
         rowShift = cRow*blockRows;
         colShift = cRow*blockCols;
      }

      // Create real main LHS matrix
      SparseMatrix work =  linearCoeff*spEq->linearMatrix(comp,idx).first - timeCoeff*spEq->timeMatrix(comp,idx).first;

      // reserve storage for the triplets
      int nz = work.nonZeros();
      nz += spEq->bcMatrix(comp,idx).first.nonZeros();
      if(nC > 0)
      {
         nz += spEq->couplingMatrix(comp,idx).first.nonZeros();
         nz += spEq->cbcMatrix(comp,idx).first.nonZeros();
      }
      triplets.reserve(nz);

      // Add triplets for the main LHS Matrix
      this->addTriplets(triplets, work, rowShift, colShift, 1.0);

      // Add triplets for the real boundary conditions
      this->addTriplets(triplets, spEq->bcMatrix(comp,idx).first, rowShift, colShift, 1.0);

      // Add coupling if required
      if(nC > 0)
      {
         // Set block shift for coupling LHS matrix
         rowShift = cRow*blockRows;
         colShift = (1-cRow)*blockCols;

         // Add triplets for the real coupling matrix
         this->addTriplets(triplets, spEq->couplingMatrix(comp,idx).first, rowShift, colShift, linearCoeff);

         // Add triplets for the real coupled boundary conditions
         this->addTriplets(triplets, spEq->cbcMatrix(comp,idx).first, rowShift, colShift, 1.0);
      }

      return  (nC+1)*blockRows;
   }

   int  Timestepper::buildLHSMatrix(std::vector<TripletZ>& triplets, Equations::SharedIEvolutionEquation spEq, FieldComponents::Spectral::Id comp, const int idx, const int nC, const int cRow)
   {
      // Matrix coefficients
      MHDFloat timeCoeff = ImExRK3::lhsT(this->mStep)*1.0/this->mDt;
      MHDFloat linearCoeff = ImExRK3::lhsL(this->mStep);

      // Set block size
      int blockRows = std::max(spEq->linearMatrix(comp,idx).first.rows(),spEq->linearMatrix(comp,idx).second.rows());
      int blockCols = std::max(spEq->linearMatrix(comp,idx).first.cols(),spEq->linearMatrix(comp,idx).second.cols());

      // Block shifts
      int rowShift = 0;
      int colShift = 0;

      // Set block shift for main LHS matrix
      if(nC > 0)
      {
         rowShift = cRow*blockRows;
         colShift = cRow*blockCols;
      }

      // Create main LHS matrix
      SparseMatrixZ work =  (linearCoeff*spEq->linearMatrix(comp,idx).first - timeCoeff*spEq->timeMatrix(comp,idx).first).cast<MHDComplex>() + MathConstants::cI*(linearCoeff*spEq->linearMatrix(comp,idx).second - timeCoeff*spEq->timeMatrix(comp,idx).second);

      // reserve storage for the triplets
      int nz = work.nonZeros();
      nz += spEq->bcMatrix(comp,idx).first.nonZeros() + spEq->bcMatrix(comp,idx).second.nonZeros();
      if(nC > 0)
      {
         nz += spEq->couplingMatrix(comp,idx).first.nonZeros() + spEq->couplingMatrix(comp,idx).second.nonZeros();
         nz += spEq->cbcMatrix(comp,idx).first.nonZeros() + spEq->cbcMatrix(comp,idx).second.nonZeros();
      }
      triplets.reserve(nz);

      // Add triplets for the main LHS Matrix
      this->addTriplets(triplets, work, rowShift, colShift, 1.0);

      // Create complex boundary condition matrix
      work = spEq->bcMatrix(comp,idx).first.cast<MHDComplex>() + MathConstants::cI*spEq->bcMatrix(comp,idx).second;

      // Add triplets for the boundary conditions
      this->addTriplets(triplets, work, rowShift, colShift, 1.0);

      // Add coupling if required
      if(nC > 0)
      {
         // Set block shift for coupling LHS matrix
         rowShift = cRow*blockRows;
         colShift = (1-cRow)*blockCols;

         // Create complex coupling matrix
         work = spEq->couplingMatrix(comp,idx).first.cast<MHDComplex>() + MathConstants::cI*spEq->couplingMatrix(comp,idx).second;

         // Add triplets for the real coupling matrix
         this->addTriplets(triplets, work, rowShift, colShift, linearCoeff);

         // Create complex coupled boundary condition matrix
         work = spEq->cbcMatrix(comp,idx).first.cast<MHDComplex>() + MathConstants::cI*spEq->cbcMatrix(comp,idx).second;

         // Add triplets for the real coupled boundary conditions
         this->addTriplets(triplets, work, rowShift, colShift, 1.0);
      }

      return  (nC+1)*blockRows;
   }

   void  Timestepper::buildRHSMatrix(std::vector<Triplet>& triplets, Equations::SharedIEvolutionEquation spEq, FieldComponents::Spectral::Id comp, const int idx, const int nC, const int cRow)
   {
      // Matrix coefficients
      MHDFloat timeCoeff = -ImExRK3::rhsT(this->mStep)*1.0/this->mDt;
      MHDFloat linearCoeff = -ImExRK3::rhsL(this->mStep);

      // Set block size
      int blockRows = spEq->linearMatrix(comp,idx).first.rows();
      int blockCols = spEq->linearMatrix(comp,idx).first.cols();

      // Block shifts
      int rowShift = 0;
      int colShift = 0;

      // Set block shift for main RHS matrix
      if(nC > 0)
      {
         rowShift = cRow*blockRows;
         colShift = cRow*blockCols;
      }

      // Create real main RHS matrix
      SparseMatrix work =  linearCoeff*spEq->linearMatrix(comp,idx).first + timeCoeff*spEq->timeMatrix(comp,idx).first;

      // reserve storage for the triplets
      int nz = work.nonZeros();
      if(nC > 0)
      {
         nz += spEq->couplingMatrix(comp,idx).first.nonZeros();
      }
      triplets.reserve(nz);

      // Add triplets for the main RHS Matrix
      this->addTriplets(triplets, work, rowShift, colShift, 1.0);

      // Add coupling if required
      if(nC > 0)
      {
         // Set block shift for coupling RHS matrix
         rowShift = cRow*blockRows;
         colShift = (1-cRow)*blockCols;

         // Add triplets for the real coupling matrix
         this->addTriplets(triplets, spEq->couplingMatrix(comp,idx).first, rowShift, colShift, linearCoeff);
      }
   }

   void  Timestepper::buildRHSMatrix(std::vector<TripletZ>& triplets, Equations::SharedIEvolutionEquation spEq, FieldComponents::Spectral::Id comp, const int idx, const int nC, const int cRow)
   {
      // Matrix coefficients
      MHDFloat timeCoeff = -ImExRK3::rhsT(this->mStep)*1.0/this->mDt;
      MHDFloat linearCoeff = -ImExRK3::rhsL(this->mStep);

      // Set block size
      int blockRows = std::max(spEq->linearMatrix(comp,idx).first.rows(),spEq->linearMatrix(comp,idx).second.rows());
      int blockCols = std::max(spEq->linearMatrix(comp,idx).first.cols(),spEq->linearMatrix(comp,idx).second.cols());

      // Block shifts
      int rowShift = 0;
      int colShift = 0;

      // Set block shift for main RHS matrix
      if(nC > 0)
      {
         rowShift = cRow*blockRows;
         colShift = cRow*blockCols;
      }

      // Create main RHS matrix
      SparseMatrixZ work =  (linearCoeff*spEq->linearMatrix(comp,idx).first + timeCoeff*spEq->timeMatrix(comp,idx).first).cast<MHDComplex>() + MathConstants::cI*(linearCoeff*spEq->linearMatrix(comp,idx).second + timeCoeff*spEq->timeMatrix(comp,idx).second);

      // reserve storage for the triplets
      int nz = work.nonZeros();
      if(nC > 0)
      {
         nz += spEq->couplingMatrix(comp,idx).first.nonZeros() + spEq->couplingMatrix(comp,idx).second.nonZeros();
      }
      triplets.reserve(nz);

      // Add triplets for the main RHS Matrix
      this->addTriplets(triplets, work, rowShift, colShift, 1.0);

      // Add coupling if required
      if(nC > 0)
      {
         // Set block shift for coupling RHS matrix
         rowShift = cRow*blockRows;
         colShift = (1-cRow)*blockCols;

         // Create complex coupling matrix
         work = spEq->couplingMatrix(comp,idx).first.cast<MHDComplex>() + MathConstants::cI*spEq->couplingMatrix(comp,idx).second;

         // Add triplets for the coupling matrix
         this->addTriplets(triplets, work, rowShift, colShift, linearCoeff);
      }
   }

   void Timestepper::addTriplets(std::vector<Triplet>& triplets, const SparseMatrix& mat, const int rowShift, const int colShift, const MHDFloat c)
   {
      // Add triplet for matrix
      for (int k=0; k<mat.outerSize(); ++k)
      {
         for (SparseMatrix::InnerIterator it(mat,k); it; ++it)
         {
            triplets.push_back(Triplet(it.row()+rowShift, it.col()+colShift, c*it.value()));
         }
      }
   }

   void Timestepper::addTriplets(std::vector<TripletZ>& triplets, const SparseMatrixZ& mat, const int rowShift, const int colShift, const MHDComplex c)
   {
      // Add triplet for matrix
      for (int k=0; k<mat.outerSize(); ++k)
      {
         for (SparseMatrixZ::InnerIterator it(mat,k); it; ++it)
         {
            triplets.push_back(TripletZ(it.row()+rowShift, it.col()+colShift, c*it.value()));
         }
      }
   }

   void  Timestepper::updateLHSMatrix(std::vector<Triplet>& triplets, Equations::SharedIEvolutionEquation spEq, FieldComponents::Spectral::Id comp, const int idx, const int nC, const int cRow)
   {
      // Matrix coefficients
      MHDFloat timeCoeff = ImExRK3::lhsT(this->mStep)*(1.0/this->mOldDt - 1.0/this->mDt);

      // Set block size
      int blockRows = spEq->timeMatrix(comp,idx).first.rows();
      int blockCols = spEq->timeMatrix(comp,idx).first.cols();

      // Block shifts
      int rowShift = 0;
      int colShift = 0;

      // Set block shift for time matrix
      if(nC > 0)
      {
         rowShift = cRow*blockRows;
         colShift = cRow*blockCols;
      }

      // Create main LHS matrix
      SparseMatrix work =  timeCoeff*spEq->timeMatrix(comp,idx).first;

      std::vector<Triplet> timeTriplets;
      timeTriplets.reserve(work.nonZeros());

      // Add triplets for the main LHS Matrix
      this->addTriplets(timeTriplets, work, rowShift, colShift, 1.0);

      std::vector<Triplet>::iterator it;
      size_t j = 0;
      for(it = triplets.begin(); it != triplets.end(); ++it)
      {
         if(it->col() == timeTriplets.at(j).col() && it->row() == timeTriplets.at(j).row())
         {
            *it = Triplet(it->row(), it->col(), it->value() + timeTriplets.at(j).value());
            j++;

            if(j == timeTriplets.size())
            {
               break;
            }
         }
      }

      // Safety assert that all values have been updated
      assert(j == timeTriplets.size());
   }

   void  Timestepper::updateLHSMatrix(std::vector<TripletZ>& triplets, Equations::SharedIEvolutionEquation spEq, FieldComponents::Spectral::Id comp, const int idx, const int nC, const int cRow)
   {
      // Matrix coefficients
      MHDFloat timeCoeff = ImExRK3::lhsT(this->mStep)*(1.0/this->mOldDt - 1.0/this->mDt);

      // Set block size
      int blockRows = std::max(spEq->timeMatrix(comp,idx).first.rows(),spEq->timeMatrix(comp,idx).second.rows());
      int blockCols = std::max(spEq->timeMatrix(comp,idx).first.cols(),spEq->timeMatrix(comp,idx).second.cols());

      // Block shifts
      int rowShift = 0;
      int colShift = 0;

      // Set block shift for time matrix
      if(nC > 0)
      {
         rowShift = cRow*blockRows;
         colShift = cRow*blockCols;
      }

      // Create main LHS matrix
      SparseMatrixZ work =  timeCoeff*spEq->timeMatrix(comp,idx).first.cast<MHDComplex>() + MathConstants::cI*timeCoeff*spEq->timeMatrix(comp,idx).second;

      std::vector<TripletZ> timeTriplets;
      timeTriplets.reserve(work.nonZeros());

      // Add triplets for the main LHS Matrix
      this->addTriplets(timeTriplets, work, rowShift, colShift, 1.0);

      std::vector<TripletZ>::iterator it;
      size_t j = 0;
      for(it = triplets.begin(); it != triplets.end(); ++it)
      {
         if(it->col() == timeTriplets.at(j).col() && it->row() == timeTriplets.at(j).row())
         {
            *it = TripletZ(it->row(), it->col(), it->value() + timeTriplets.at(j).value());
            j++;

            if(j == timeTriplets.size())
            {
               break;
            }
         }
      }

      // Safety assert that all values have been updated
      assert(j == timeTriplets.size());
   }

   void  Timestepper::updateRHSMatrix(std::vector<Triplet>& triplets, Equations::SharedIEvolutionEquation spEq, FieldComponents::Spectral::Id comp, const int idx, const int nC, const int cRow)
   {
      // Matrix coefficients
      MHDFloat timeCoeff = ImExRK3::lhsT(this->mStep)*(1.0/this->mOldDt - 1.0/this->mDt);

      // Set block size
      int blockRows = spEq->timeMatrix(comp,idx).first.rows();
      int blockCols = spEq->timeMatrix(comp,idx).first.cols();

      // Block shifts
      int rowShift = 0;
      int colShift = 0;

      // Set block shift for time matrix
      if(nC > 0)
      {
         rowShift = cRow*blockRows;
         colShift = cRow*blockCols;
      }

      // Create main LHS matrix
      SparseMatrix work =  timeCoeff*spEq->timeMatrix(comp,idx).first;

      std::vector<Triplet> timeTriplets;
      timeTriplets.reserve(work.nonZeros());

      // Add triplets for the main LHS Matrix
      this->addTriplets(timeTriplets, work, rowShift, colShift, 1.0);

      std::vector<Triplet>::iterator it;
      size_t j = 0;
      for(it = triplets.begin(); it != triplets.end(); ++it)
      {
         if(it->col() == timeTriplets.at(j).col() && it->row() == timeTriplets.at(j).row())
         {
            *it = Triplet(it->row(), it->col(), it->value() + timeTriplets.at(j).value());
            j++;

            if(j == timeTriplets.size())
            {
               break;
            }
         }
      }

      // Safety assert that all values have been updated
      assert(j == timeTriplets.size());
   }

   void  Timestepper::updateRHSMatrix(std::vector<TripletZ>& triplets, Equations::SharedIEvolutionEquation spEq, FieldComponents::Spectral::Id comp, const int idx, const int nC, const int cRow)
   {
      // Matrix coefficients
      MHDFloat timeCoeff = ImExRK3::rhsT(this->mStep)*(1.0/this->mOldDt - 1.0/this->mDt);

      // Set block size
      int blockRows = std::max(spEq->timeMatrix(comp,idx).first.rows(),spEq->timeMatrix(comp,idx).second.rows());
      int blockCols = std::max(spEq->timeMatrix(comp,idx).first.cols(),spEq->timeMatrix(comp,idx).second.cols());

      // Block shifts
      int rowShift = 0;
      int colShift = 0;

      // Set block shift for time matrix
      if(nC > 0)
      {
         rowShift = cRow*blockRows;
         colShift = cRow*blockCols;
      }

      // Create main LHS matrix
      SparseMatrixZ work =  timeCoeff*spEq->timeMatrix(comp,idx).first.cast<MHDComplex>() + MathConstants::cI*timeCoeff*spEq->timeMatrix(comp,idx).second;

      std::vector<TripletZ> timeTriplets;
      timeTriplets.reserve(work.nonZeros());

      // Add triplets for the main LHS Matrix
      this->addTriplets(timeTriplets, work, rowShift, colShift, 1.0);

      std::vector<TripletZ>::iterator it;
      size_t j = 0;
      for(it = triplets.begin(); it != triplets.end(); ++it)
      {
         if(it->col() == timeTriplets.at(j).col() && it->row() == timeTriplets.at(j).row())
         {
            *it = TripletZ(it->row(), it->col(), it->value() + timeTriplets.at(j).value());
            j++;

            if(j == timeTriplets.size())
            {
               break;
            }
         }
      }

      // Safety assert that all values have been updated
      assert(j == timeTriplets.size());
   }
}
}
