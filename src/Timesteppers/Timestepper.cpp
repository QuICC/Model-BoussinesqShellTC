/** \file Timestepper.cpp
 *  \brief Implementation of a general timestepper structure
 */

// Debug includes
//
#include "Debug/DebuggerMacro.h"

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
         DebuggerMacro_showValue("Updating timestep and matrices with new Dt = ", 0, this->mDt);
      }

      //
      // Update the timestep matrices if necessary
      //
      if(hasNewDt)
      {
         DebuggerMacro_start("Update matrices", 0);
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
         DebuggerMacro_stop("Update matrices t = ", 0);

         DebuggerMacro_start("Complex solver update", 0);
         // Update solvers from complex equation steppers
         for(size_t i = 0; i < this->mEqZStepper.size(); i++)
         {
            this->mEqZStepper.at(i).updateSolver();
         }
         DebuggerMacro_stop("Complex solver update t = ", 0);

         DebuggerMacro_start("Real solver update", 0);
         // Update solvers from real equation steppers
         for(size_t i = 0; i < this->mEqDStepper.size(); i++)
         {
            this->mEqDStepper.at(i).updateSolver();
         }
         DebuggerMacro_stop("Real solver update t = ", 0);

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

      DebuggerMacro_start("Create matrices", 0);
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
      DebuggerMacro_stop("Create matrices t = ", 0);

      //
      // Initialise the solvers and the initial state
      //

      DebuggerMacro_start("Complex solver init", 0);
      // Initialise solvers from complex equation steppers
      for(size_t i = 0; i < this->mEqZStepper.size(); i++)
      {
         this->mEqZStepper.at(i).initSolver();
      }
      DebuggerMacro_stop("Complex solver init t = ", 0);

      DebuggerMacro_start("Real solver init", 0);
      // Initialise solvers from real equation steppers
      for(size_t i = 0; i < this->mEqDStepper.size(); i++)
      {
         this->mEqDStepper.at(i).initSolver();
      }
      DebuggerMacro_stop("Real solver init t = ", 0);

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
               // Add LHS triplets
               size = this->buildSolverMatrix(this->mEqZStepper.at(myIdx).rLHSMatrix(start+i), spEq, comp, i, nC, fieldIdx, true);

               // Add RHS triplets
               this->buildSolverMatrix(this->mEqZStepper.at(myIdx).rRHSMatrix(start+i), spEq, comp, i, nC, fieldIdx, false);

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
               this->buildSolverMatrix(this->mEqZStepper.at(myIdx).rLHSMatrix(start+i), spEq, comp, i, nC, fieldIdx, true);

               // Add RHS triplets
               this->buildSolverMatrix(this->mEqZStepper.at(myIdx).rRHSMatrix(start+i), spEq, comp, i, nC, fieldIdx, false);

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
               // Add LHS triplets
               size = this->buildSolverMatrix(this->mEqDStepper.at(myIdx).rLHSMatrix(start+i), spEq, comp, i, nC, fieldIdx, true);

               // Add RHS triplets
               this->buildSolverMatrix(this->mEqDStepper.at(myIdx).rRHSMatrix(start+i), spEq, comp, i, nC, fieldIdx, false);

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
               this->buildSolverMatrix(this->mEqDStepper.at(myIdx).rLHSMatrix(start+i), spEq, comp, i, nC, fieldIdx, true);

               // Add RHS triplets
               this->buildSolverMatrix(this->mEqDStepper.at(myIdx).rRHSMatrix(start+i), spEq, comp, i, nC, fieldIdx, false);

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
               this->updateTimeMatrix(this->mEqZStepper.at(myIdx).rLHSMatrix(start+i), spEq, comp, i, nC, fieldIdx, true);

               // Update RHS triplets
               this->updateTimeMatrix(this->mEqZStepper.at(myIdx).rRHSMatrix(start+i), spEq, comp, i, nC, fieldIdx, false);

               // Set the start row
               startRow(i) = 0;
            }
         } else
         {
            for(int i = 0; i < interInfo.first; i++)
            {
               // Update LHS triplets
               this->updateTimeMatrix(this->mEqZStepper.at(myIdx).rLHSMatrix(start+i), spEq, comp, i, nC, fieldIdx, true);

               // Update RHS triplets
               this->updateTimeMatrix(this->mEqZStepper.at(myIdx).rRHSMatrix(start+i), spEq, comp, i, nC, fieldIdx, false);

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
               this->updateTimeMatrix(this->mEqDStepper.at(myIdx).rLHSMatrix(start+i), spEq, comp, i, nC, fieldIdx, true);

               // Update RHS triplets
               this->updateTimeMatrix(this->mEqDStepper.at(myIdx).rRHSMatrix(start+i), spEq, comp, i, nC, fieldIdx, false);

               // Set the start row
               startRow(i) = 0;
            }
         } else
         {
            for(int i = 0; i < interInfo.first; i++)
            {
               // Update LHS triplets
               this->updateTimeMatrix(this->mEqDStepper.at(myIdx).rLHSMatrix(start+i), spEq, comp, i, nC, fieldIdx, true);

               // Update RHS triplets
               this->updateTimeMatrix(this->mEqDStepper.at(myIdx).rRHSMatrix(start+i), spEq, comp, i, nC, fieldIdx, false);

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

   int Timestepper::buildSolverMatrix(SparseMatrix& solverMatrix, Equations::SharedIEvolutionEquation spEq, FieldComponents::Spectral::Id comp, const int idx, const int nC, const int cRow, const bool isLhs)
   {
      // Operator coefficients
      MHDFloat timeCoeff;
      MHDFloat linearCoeff;

      if(isLhs)
      {
         // Set time coefficients for LHS Matrix
         timeCoeff = ImExRK3::lhsT(this->mStep)*1.0/this->mDt;

         // Set linear coefficients for LHS Matrix
         linearCoeff = ImExRK3::lhsL(this->mStep);
      } else
      {
         // Set time coefficients for RHS Matrix
         timeCoeff = ImExRK3::rhsT(this->mStep)*1.0/this->mDt;

         // Set linear coefficients for RHS Matrix
         linearCoeff = -ImExRK3::rhsL(this->mStep);
      }

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

      // Create internal part of solver matrix (i.e without boundary conditions)
      SparseMatrix work =  linearCoeff*spEq->linearMatrix(comp,idx).first - timeCoeff*spEq->timeMatrix(comp,idx).first;

      // Count the number of nonzero values for internal part of solver matrix
      int nz = work.nonZeros();
      if(nC > 0)
      {
         nz += spEq->couplingMatrix(comp,idx).first.nonZeros();
      }

      // Count the number of nonzero values for boundary conditions
      if(isLhs)
      {
         nz += spEq->bcMatrix(comp,idx).first.nonZeros();
         if(nC > 0)
         {
            nz += spEq->cbcMatrix(comp,idx).first.nonZeros();
         }
      }
      
      // Storage for the triplets
      std::vector<Triplet> triplets;
      // Reserve space to avoid reallocations
      triplets.reserve(nz + solverMatrix.nonZeros());
      // Add nonzero elements from solver matrix to triplets
      this->addTriplets(triplets, solverMatrix, 0, 0, 1.0);

      // Add triplets for the internal part of solver Matrix
      this->addTriplets(triplets, work, rowShift, colShift, 1.0);

      // Add coupling if required
      if(nC > 0)
      {
         // Set block shift for coupling matrix
         rowShift = cRow*blockRows;
         colShift = (1-cRow)*blockCols;

         // Add triplets for the coupling matrix
         this->addTriplets(triplets, spEq->couplingMatrix(comp,idx).first, rowShift, colShift, linearCoeff);
      }

      if(isLhs)
      {
         // Set block shift
         rowShift = cRow*blockRows;
         colShift = cRow*blockCols;

         // Add triplets for the boundary conditions
         this->addTriplets(triplets, spEq->bcMatrix(comp,idx).first, rowShift, colShift, 1.0);

         // Add coupling if required
         if(nC > 0)
         {
            // Set block shift for coupling matrix
            rowShift = cRow*blockRows;
            colShift = (1-cRow)*blockCols;

            // Add triplets for the coupled boundary conditions
            this->addTriplets(triplets, spEq->cbcMatrix(comp,idx).first, rowShift, colShift, 1.0);
         }
      }

      // Set matrix from triplets
      solverMatrix.resize((nC+1)*blockRows,(nC+1)*blockRows);
      solverMatrix.setFromTriplets(triplets.begin(), triplets.end());

      return  solverMatrix.rows();
   }

   int  Timestepper::buildSolverMatrix(SparseMatrixZ& solverMatrix, Equations::SharedIEvolutionEquation spEq, FieldComponents::Spectral::Id comp, const int idx, const int nC, const int cRow, const bool isLhs)
   {
      // Operator coefficients
      MHDFloat timeCoeff;
      MHDFloat linearCoeff;

      if(isLhs)
      {
         // Set time coefficients for LHS Matrix
         timeCoeff = ImExRK3::lhsT(this->mStep)*1.0/this->mDt;

         // Set linear coefficients for LHS Matrix
         linearCoeff = ImExRK3::lhsL(this->mStep);
      } else
      {
         // Set time coefficients for RHS Matrix
         timeCoeff = ImExRK3::rhsT(this->mStep)*1.0/this->mDt;

         // Set linear coefficients for RHS Matrix
         linearCoeff = -ImExRK3::rhsL(this->mStep);
      }

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

      // Count the number of nonzero values for internal part of solver matrix
      int nz = work.nonZeros();
      if(nC > 0)
      {
         nz += spEq->couplingMatrix(comp,idx).first.nonZeros() + spEq->couplingMatrix(comp,idx).second.nonZeros();
      }

      // Count the number of nonzero values for boundary conditions
      if(isLhs)
      {
         nz += spEq->bcMatrix(comp,idx).first.nonZeros() + spEq->bcMatrix(comp,idx).second.nonZeros();
         if(nC > 0)
         {
            nz += spEq->cbcMatrix(comp,idx).first.nonZeros() + spEq->cbcMatrix(comp,idx).second.nonZeros();
         }
      }
      
      // Storage for the triplets
      std::vector<TripletZ> triplets;
      // Reserve space to avoid reallocations
      triplets.reserve(nz + solverMatrix.nonZeros());
      // Add nonzero elements from solver matrix to triplets
      this->addTriplets(triplets, solverMatrix, 0, 0, 1.0);

      // Add triplets for the main LHS Matrix
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
      }

      if(isLhs)
      {
         // Set block shift
         rowShift = cRow*blockRows;
         colShift = cRow*blockCols;

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

            // Create complex coupled boundary condition matrix
            work = spEq->cbcMatrix(comp,idx).first.cast<MHDComplex>() + MathConstants::cI*spEq->cbcMatrix(comp,idx).second;

            // Add triplets for the real coupled boundary conditions
            this->addTriplets(triplets, work, rowShift, colShift, 1.0);
         }
      }

      // Set matrix from triplets
      solverMatrix.resize((nC+1)*blockRows,(nC+1)*blockRows);
      solverMatrix.setFromTriplets(triplets.begin(), triplets.end());

      return  solverMatrix.rows();
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

   void  Timestepper::updateTimeMatrix(SparseMatrix& oldTime, Equations::SharedIEvolutionEquation spEq, FieldComponents::Spectral::Id comp, const int idx, const int nC, const int cRow, bool isLhs)
   {
      // Set time coefficient
      MHDFloat timeCoeff;
      if(isLhs)
      {
         // Compute timestep correction coefficient for LHS matrix
         timeCoeff = ImExRK3::lhsT(this->mStep)*(1.0/this->mOldDt - 1.0/this->mDt);
      } else
      {
         // Compute timestep correction coefficient for RHS matrix
         timeCoeff = ImExRK3::rhsT(this->mStep)*(1.0/this->mOldDt - 1.0/this->mDt);
      }

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

      // Convert time matrix to triplets
      std::vector<Triplet> timeTriplets;
      timeTriplets.reserve(work.nonZeros());
      this->addTriplets(timeTriplets, work, rowShift, colShift, 1.0);

      // Update old matrix with new time dependence
      int j = 0;
      for (int k=0; k<oldTime.outerSize(); ++k)
      {
         for (SparseMatrix::InnerIterator it(oldTime,k); it; ++it)
         {
            if(static_cast<size_t>(j) < timeTriplets.size() && static_cast<size_t>(k) == timeTriplets.at(j).col())
            {
               if(static_cast<size_t>(it.row()) == timeTriplets.at(j).row())
               {
                  it.valueRef() += timeTriplets.at(j).value();
                  j++;
               }
            }
         }
      }

      // Safety assert to make sure all values have been updated
      assert(static_cast<size_t>(j) == timeTriplets.size());
   }

   void  Timestepper::updateTimeMatrix(SparseMatrixZ& oldTime, Equations::SharedIEvolutionEquation spEq, FieldComponents::Spectral::Id comp, const int idx, const int nC, const int cRow, bool isLhs)
   {
      // Set time coefficient
      MHDFloat timeCoeff;
      if(isLhs)
      {
         // Compute timestep correction coefficient for LHS matrix
         timeCoeff = ImExRK3::lhsT(this->mStep)*(1.0/this->mOldDt - 1.0/this->mDt);
      } else
      {
         // Compute timestep correction coefficient for RHS matrix
         timeCoeff = ImExRK3::rhsT(this->mStep)*(1.0/this->mOldDt - 1.0/this->mDt);
      }

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

      // Create time matrix
      SparseMatrixZ work =  timeCoeff*spEq->timeMatrix(comp,idx).first.cast<MHDComplex>() + MathConstants::cI*timeCoeff*spEq->timeMatrix(comp,idx).second;

      // Convert time matrix to triplets
      std::vector<TripletZ> timeTriplets;
      timeTriplets.reserve(work.nonZeros());
      this->addTriplets(timeTriplets, work, rowShift, colShift, 1.0);

      // Update old matrix with new time dependence
      int j = 0;
      for (int k=0; k<oldTime.outerSize(); ++k)
      {
         for (SparseMatrixZ::InnerIterator it(oldTime,k); it; ++it)
         {
            if(static_cast<size_t>(j) < timeTriplets.size() && static_cast<size_t>(k) == timeTriplets.at(j).col())
            {
               if(static_cast<size_t>(it.row()) == timeTriplets.at(j).row())
               {
                  it.valueRef() += timeTriplets.at(j).value();
                  j++;
               }
            }
         }
      }

      // Safety assert to make sure all values have been updated
      assert(static_cast<size_t>(j) == timeTriplets.size());
   }
}
}
