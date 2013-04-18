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
      // Create real/complex timesteppers
      //

      // Loop over all scalar equations
      std::vector<Equations::SharedIScalarEquation>::const_iterator scalEqIt;
      for(scalEqIt = scalEq.begin(); scalEqIt < scalEq.end(); scalEqIt++)
      {
         // Get type information for the equation steppers
         this->createEqStepper((*scalEqIt), FieldComponents::Spectral::SCALAR);
      }

      // Loop over all vector equations
      std::vector<Equations::SharedIVectorEquation>::const_iterator vectEqIt;
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

   void Timestepper::createEqStepper(Equations::SharedIEvolutionEquation spEq, FieldComponents::Spectral::Id comp)
   {
      // Equation is part of a complex system
      if(spEq->couplingInfo(comp).isComplex())
      {
         // Add equation stepper if system index does not yet exist
         if(spEq->couplingInfo(comp).solverIndex() > this->mEqZStepper.size() - 1)
         {
            this->mEqZStepper.push_back(EquationZTimestepper(spEq->systemFields(comp), spEq->startIndex(comp)));
         }

      // Equation is part of a real system
      } else
      {
         // Add equation stepper if system index does not yet exist
         if(spEq->couplingInfo(comp).solverIndex() > this->mEqDStepper.size() - 1)
         {
            this->mEqDStepper.push_back(EquationZTimestepper(spEq->systemFields(comp), spEq->startIndex(comp)));
         }
      }
   }

   void Timestepper::createMatrices(Equations::SharedIEvolutionEquation spEq, FieldComponents::Spectral::Id comp)
   {
      // ID of the current field
      std::pair<PhysicalNames::Id, FieldComponents::Spectral::Id> myId = std::make_pair(spEq->name(),comp);

      // Index of the current field
      int myIdx = spEq->couplingInfo(comp).solverIndex();

      // Number of linear systems
      int nSystems = spEq->couplingInfo(comp).nSystems();

      // start index for matrices
      int start = this->mStep*nSystems;

      // Start row for storage information
      ArrayI startRow(nSystems);

      // Complex matrices in linear solve
      if(spEq->couplingInfo(comp).isComplex())
      {
        // Reserve space for the matrices to avoid large number of expensive reallocations
        if(spEq->couplingInfo(comp).solverIndex() == 0 && this->mStep == 0)
        {
           this->mEqZStepper.at(myIdx).reserveMatrices(ImExRK3::STEPS*nSystems);
        }

        int size = 0;
        for(int i = 0; i < nSystems; i++)
        {
           // Add LHS triplets
           this->buildSolverMatrix(this->mEqZStepper.at(myIdx).rLHSMatrix(start+i), spEq, comp, i, true);

           // Add RHS triplets
           this->buildSolverMatrix(this->mEqZStepper.at(myIdx).rRHSMatrix(start+i), spEq, comp, i, false);

           if(spEq->systemIndex(comp) == 0 && this->mStep == 0)
           {
              // Create RHS and solution data storage
              this->mEqZStepper.at(myIdx).addStorage(size, spEq->couplingInfo(comp).rhsCols(i));
           }

           // Set the start row
           startRow(i) = spEq->systemIndex(comp)*size;
        }

         if(this->mStep == 0)
         {
            // Store storage information
            this->mEqZStepper.at(myIdx).addInformation(myId,startRow);
         }

      // Real matrices in linear solve
      } else
      {
         // Reserve space for the matrices to avoid large number of expensive reallocations
         if(spEq->systemIndex(comp) == 0 && this->mStep == 0)
         {
            this->mEqDStepper.at(myIdx).reserveMatrices(ImExRK3::STEPS*nSystems);
         }

         int size = 0;
         for(int i = 0; i < nSystems.first; i++)
         {
            // Add LHS triplets
            this->buildSolverMatrix(this->mEqDStepper.at(myIdx).rLHSMatrix(start+i), spEq, comp, i, true);

            // Add RHS triplets
            this->buildSolverMatrix(this->mEqDStepper.at(myIdx).rRHSMatrix(start+i), spEq, comp, i, false);

            if(spEq->systemIndex(comp) == 0 && this->mStep == 0)
            {
               // Create RHS and solution data storage
               this->mEqDStepper.at(myIdx).addStorage(size, spEq->couplingInfo(comp).rhsCols(i));
            }

            // Set the start row
            startRow(i) = spEq->systemIndex(comp)*size;
         }

         if(this->mStep == 0)
         {
            // Store storage information
            this->mEqDStepper.at(myIdx).addInformation(myId,startRow);
         }
      }
   }

   void Timestepper::updateMatrices(Equations::SharedIEvolutionEquation spEq, FieldComponents::Spectral::Id comp)
   {
      // ID of the current field
      std::pair<PhysicalNames::Id, FieldComponents::Spectral::Id> myId = std::make_pair(spEq->name(),comp);

      // Index of the current field
      int myIdx = spEq->couplingInfo(comp).solverIndex();

      // start index for matrices
      int start = this->mStep*interInfo.first;

      // Complex matrices in linear solve
      if(spEq->couplingInfo(comp).isComplex())
      {
         for(int i = 0; i < interInfo.first; i++)
         {
            // Update LHS triplets
            this->updateTimeMatrix(this->mEqZStepper.at(myIdx).rLHSMatrix(start+i), spEq, comp, i, true);

            // Update RHS triplets
            this->updateTimeMatrix(this->mEqZStepper.at(myIdx).rRHSMatrix(start+i), spEq, comp, i, false);
         }

      // Real matrices in linear solve
      } else
      {
         for(int i = 0; i < interInfo.first; i++)
         {
            // Update LHS triplets
            this->updateTimeMatrix(this->mEqDStepper.at(myIdx).rLHSMatrix(start+i), spEq, comp, i, true);

            // Update RHS triplets
            this->updateTimeMatrix(this->mEqDStepper.at(myIdx).rRHSMatrix(start+i), spEq, comp, i, false);
         }
      }
   }

   void Timestepper::initSolution(const std::vector<Equations::SharedIScalarEquation>& scalEq, const std::vector<Equations::SharedIVectorEquation>& vectEq)
   {
      // Storage for information and identity
      std::pair<PhysicalNames::Id, FieldComponents::Spectral::Id> myId;

      // Storage for the selected field component
      FieldComponents::Spectral::Id comp;

      // Loop over all scalar equations
      comp = FieldComponents::Spectral::SCALAR;
      std::vector<Equations::SharedIScalarEquation>::const_iterator scalEqIt;
      for(scalEqIt = scalEq.begin(); scalEqIt < scalEq.end(); scalEqIt++)
      {
         // Get identity and corresponding equation information
         myId = std::make_pair((*scalEqIt)->name(), comp);

         // Get index of current field
         int myIdx = (*scalEqIt)->couplingInfo(comp).solverIndex();

         // Linear solve matrices are complex
         if((*scalEqIt)->couplingInfo(comp).isComplex())
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
         // Work on first component
         comp = FieldComponents::Spectral::ONE;

         // Get identity and corresponding equation information for toroidal component
         myId = std::make_pair((*vectEqIt)->name(), comp);

         // Get index of current field
         int myIdx = (*vectEqIt)->couplingInfo(comp).solverIndex();

         // Linear solve matrices are complex
         if((*vectEqIt)->couplingInfo(comp).isComplex())
         {
            // Get timestep input for toroidal component
            for(int i = 0; i < this->mEqZStepper.at(myIdx).nSystem(); i++)
            {
               (*vectEqIt)->copyTInput(comp, this->mEqZStepper.at(myIdx).rSolution(i), i, this->mEqZStepper.at(myIdx).startRow(myId,i));
            }

         // Linear solve matrices are real
         } else
         {
            // Get timestep input for toroidal component
            for(int i = 0; i < this->mEqDStepper.at(myIdx).nSystem(); i++)
            {
               (*vectEqIt)->copyTInput(comp, this->mEqDStepper.at(myIdx).rSolution(i), i, this->mEqDStepper.at(myIdx).startRow(myId,i));
            }
         }

         // Work on second component
         comp = FieldComponents::Spectral::TWO;

         // Get identity and corresponding equation information for second component
         myId = std::make_pair((*vectEqIt)->name(), comp);

         // Get index of current field
         myIdx = (*vectEqIt)->couplingInfo(comp).solverIndex();

         // Linear solve matrices are complex
         if((*vectEqIt)->couplingInfo(comp).isComplex())
         {
            // Get timestep input for poloidal component
            for(int i = 0; i < this->mEqZStepper.at(myIdx).nSystem(); i++)
            {
               (*vectEqIt)->copyTInput(comp, this->mEqZStepper.at(myIdx).rSolution(i), i, this->mEqZStepper.at(myIdx).startRow(myId,i));
            }

         // Linear solve matrices are real
         } else
         {
            // Get timestep input for poloidal component
            for(int i = 0; i < this->mEqDStepper.at(myIdx).nSystem(); i++)
            {
               (*vectEqIt)->copyTInput(comp, this->mEqDStepper.at(myIdx).rSolution(i), i, this->mEqDStepper.at(myIdx).startRow(myId,i));
            }
         }
      }
   }

   void Timestepper::getInput(const std::vector<Equations::SharedIScalarEquation>& scalEq, const std::vector<Equations::SharedIVectorEquation>& vectEq)
   {
      // Storage for information and identity
      std::pair<PhysicalNames::Id, FieldComponents::Spectral::Id> myId;

      // Storage for the selected field component
      FieldComponents::Spectral::Id comp;

      // Loop over all scalar equations
      comp = FieldComponents::Spectral::SCALAR;
      std::vector<Equations::SharedIScalarEquation>::const_iterator scalEqIt;
      for(scalEqIt = scalEq.begin(); scalEqIt < scalEq.end(); scalEqIt++)
      {
         // Get identity and corresponding equation information
         myId = std::make_pair((*scalEqIt)->name(), comp);

         // Get index of current field
         int myIdx = (*scalEqIt)->systemIndex(comp);

         // Linear solve matrices are complex
         if((*scalEqIt)->isSystemComplex(comp))
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
         comp = FieldComponents::Spectral::ONE;

         // Get identity and corresponding equation information for toroidal component
         myId = std::make_pair((*vectEqIt)->name(), comp);

         // Get index of current field
         int myIdx = (*vectEqIt)->systemIndex(comp);

         // Linear solve matrices are complex
         if((*vectEqIt)->isSystemComplex(comp))
         {
            // Get timestep input for toroidal component
            for(int i = 0; i < this->mEqZStepper.at(myIdx).nSystem(); i++)
            {
               (*vectEqIt)->timestepInput(comp, this->mEqZStepper.at(myIdx).rRHSData(i), i, this->mEqZStepper.at(myIdx).startRow(myId,i));
            }

         // Linear solve matrices are real
         } else
         {
            // Get timestep input for toroidal component
            for(int i = 0; i < this->mEqDStepper.at(myIdx).nSystem(); i++)
            {
               (*vectEqIt)->timestepInput(comp, this->mEqDStepper.at(myIdx).rRHSData(i), i, this->mEqDStepper.at(myIdx).startRow(myId,i));
            }
         }

         comp = FieldComponents::Spectral::TWO;

         // Get identity and corresponding equation information for poloidal component
         myId = std::make_pair((*vectEqIt)->name(), comp);

         // Get index of current field
         myIdx = (*vectEqIt)->systemIndex(comp);

         // Linear solve matrices are complex
         if(this->(*vectEqIt)->isSystemComplex(comp))
         {
            // Get timestep input for poloidal component
            for(int i = 0; i < this->mEqZStepper.at(myIdx).nSystem(); i++)
            {
               (*vectEqIt)->timestepInput(comp, this->mEqZStepper.at(myIdx).rRHSData(i), i, this->mEqZStepper.at(myIdx).startRow(myId,i));
            }

         // Linear solve matrices are real
         } else
         {
            // Get timestep input for poloidal component
            for(int i = 0; i < this->mEqDStepper.at(myIdx).nSystem(); i++)
            {
               (*vectEqIt)->timestepInput(comp, this->mEqDStepper.at(myIdx).rRHSData(i), i, this->mEqDStepper.at(myIdx).startRow(myId,i));
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
      std::pair<PhysicalNames::Id, FieldComponents::Spectral::Id> myId;

      FieldComponents::Spectral::Id comp;

      // Loop over all scalar equations
      comp = FieldComponents::Spectral::SCALAR;
      std::vector<Equations::SharedIScalarEquation>::const_iterator scalEqIt;
      for(scalEqIt = scalEq.begin(); scalEqIt < scalEq.end(); scalEqIt++)
      {
         // Get identity and corresponding equation information
         myId = std::make_pair((*scalEqIt)->name(), comp);

         // Get index of current field
         int myIdx = (*scalEqIt)->systemIndex(comp);

         // Linear solve matrices are complex
         if(this->(*scalEqIt)->isSystemComplex(comp))
         {
            // Get timestep output
            for(int i = 0; i < this->mEqZStepper.at(myIdx).nSystem(); i++)
            {
               (*scalEqIt)->timestepOutput(comp, this->mEqZStepper.at(myIdx).solution(i), i, this->mEqZStepper.at(myIdx).startRow(myId,i));
            }

         // Linear solve matrices are real
         } else
         {
            // Get timestep output
            for(int i = 0; i < this->mEqDStepper.at(myIdx).nSystem(); i++)
            {
               (*scalEqIt)->timestepOutput(comp, this->mEqDStepper.at(myIdx).solution(i), i, this->mEqDStepper.at(myIdx).startRow(myId,i));
            }
         }
      }

      // Loop over all vector equations
      std::vector<Equations::SharedIVectorEquation>::const_iterator vectEqIt;
      for(vectEqIt = vectEq.begin(); vectEqIt < vectEq.end(); vectEqIt++)
      {
         comp = FieldComponents::Spectral::ONE;

         // Get identity and corresponding equation information for toroidal component
         myId = std::make_pair((*vectEqIt)->name(), comp);

         // Get index of current field
         int myIdx = (*vectEqIt)->systemIndex(comp);

         // Linear solve matrices are complex
         if((*vectEqIt)->isSystemComplex(comp))
         {
            // Get timestep output for toroidal component
            for(int i = 0; i < this->mEqZStepper.at(myIdx).nSystem(); i++)
            {
               (*vectEqIt)->timestepOutput(comp, this->mEqZStepper.at(myIdx).solution(i), i, this->mEqZStepper.at(myIdx).startRow(myId,i));
            }

         // Linear solve matrices are real
         } else
         {
            // Get timestep output for toroidal component
            for(int i = 0; i < this->mEqDStepper.at(myIdx).nSystem(); i++)
            {
               (*vectEqIt)->timestepOutput(comp, this->mEqDStepper.at(myIdx).solution(i), i, this->mEqDStepper.at(myIdx).startRow(myId,i));
            }
         }

         comp = FieldComponents::Spectral::TWO;

         // Get identity and corresponding equation information for poloidal component
         myId = std::make_pair((*vectEqIt)->name(), comp);

         // Get index of current field
         myIdx = (*vectEqIt)->systemIndex(comp);

         // Linear solve matrices are complex
         if((*vectEqIt)->isSystemComplex(comp))
         {
            // Get timestep input for poloidal component
            for(int i = 0; i < this->mEqZStepper.at(myIdx).nSystem(); i++)
            {
               (*vectEqIt)->timestepOutput(comp, this->mEqZStepper.at(myIdx).solution(i), i, this->mEqZStepper.at(myIdx).startRow(myId,i));
            }

         // Linear solve matrices are real
         } else
         {
            // Get timestep output for poloidal component
            for(int i = 0; i < this->mEqDStepper.at(myIdx).nSystem(); i++)
            {
               (*vectEqIt)->timestepOutput(comp, this->mEqDStepper.at(myIdx).solution(i), i, this->mEqDStepper.at(myIdx).startRow(myId,i));
            }
         }
      }
   }

   void Timestepper::buildSolverMatrix(SparseMatrix& solverMatrix, Equations::SharedIEvolutionEquation spEq, FieldComponents::Spectral::Id comp, const int idx, const bool isLhs)
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

      // Count the number of nonzero values for linear row block
      int nz = spEq->linearRow(comp,idx).first.size();
      nz += spEq->timeRow(comp,idx).first.size();

      // Count the number of nonzero values for boundary row block
      if(isLhs)
      {
         nz += spEq->boundaryRow(comp,idx).first.size();
      }

      // Storage for the triplets
      std::vector<Triplet> solverTriplets;

      // Reserve space to avoid reallocations
      solverTriplets.reserve(nz + solverMatrix.nonZeros());

      // Add nonzero elements from solver matrix to triplets
      this->addTriplets(solverTriplets, solverMatrix, 1.0);

      // Add nonzero elements from equation linear row to triplets
      this->addTriplets(solverTriplets, spEq->linearRow(comp,idx), linearCoeff);

      // Add nonzero elements from equation linear row to triplets
      this->addTriplets(solverTriplets, spEq->timeRow(comp,idx), -timeCoeff);

      // Add boundary row for LHS operator
      if(isLhs)
      {
         this->addTriplets(solverTriplets, spEq->boundarRow(comp,idx), 1.0);
      }

      // Set matrix from triplets
      solverMatrix.resize(spEq->systemSize(comp,i),spEq->systemSize(comp,i));
      solverMatrix.setFromTriplets(solverTriplets.begin(), solverTriplets.end());
   }

   void  Timestepper::buildSolverMatrix(SparseMatrixZ& solverMatrix, Equations::SharedIEvolutionEquation spEq, FieldComponents::Spectral::Id comp, const int idx, const bool isLhs)
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

      // Count the number of nonzero values for linear row block
      int nz = spEq->linearRow(comp,idx).first.size() + spEq->linearRow(comp,idx).second.size();
      nz += spEq->timeRow(comp,idx).first.size() + spEq->timeRow(comp,idx).second.size();

      // Count the number of nonzero values for boundary row block
      if(isLhs)
      {
         nz += spEq->boundaryRow(comp,idx).first.size() + spEq->boundaryRow(comp,idx).second.size();
      }

      // Storage for the triplets
      std::vector<TripletZ> solverTriplets;

      // Reserve space to avoid reallocations
      solverTriplets.reserve(nz + solverMatrix.nonZeros());

      // Add nonzero elements from solver matrix to triplets
      this->addTriplets(solverTriplets, solverMatrix, 1.0);

      // Add nonzero elements from equation linear row to triplets
      this->addTriplets(solverTriplets, spEq->linearRow(comp,idx), linearCoeff);

      // Add nonzero elements from equation linear row to triplets
      this->addTriplets(solverTriplets, spEq->timeRow(comp,idx), -timeCoeff);

      // Add boundary row for LHS operator
      if(isLhs)
      {
         this->addTriplets(solverTriplets, spEq->boundarRow(comp,idx), 1.0);
      }

      // Set matrix from triplets
      solverMatrix.resize(spEq->systemSize(comp,i),spEq->systemSize(comp,i));
      solverMatrix.setFromTriplets(solverTriplets.begin(), solverTriplets.end());
   }

   void Timestepper::addTriplets(std::vector<Triplet>& outTriplets, const std::pair<std::vector<Triplet>,std::vector<Triplet> >& inTriplets, const MHDFloat c)
   {
      // Add real triplets to outTriplets
      for(std::vector<Triplet>::const_iterator it = inTriplets.first.begin(); it != inTriplets.first.end(); ++it)
      {
         outTriplets.push_back(Triplet(it->row(), it->col(), c*it->value()));
      }
   }

   void Timestepper::addTriplets(std::vector<TripletZ>& outTriplets, const std::pair<std::vector<Triplet>,std::vector<Triplet> >& inTriplets, const MHDFloat c)
   {
      // Add real triplets to outTriplets
      for(std::vector<Triplet>::const_iterator it = inTriplets.first.begin(); it != inTriplets.first.end(); ++it)
      {
         outTriplets.push_back(TripletZ(it->row(), it->col(), static_cast<MHDComplex>(c*it->value())));
      }

      // Add imaginary triplets to outTriplets
      for(std::vector<Triplet>::const_iterator it = inTriplets.second.begin(); it != inTriplets.second.end(); ++it)
      {
         outTriplets.push_back(TripletZ(it->row(), it->col(), c*it->value()*MathConstants::cI));
      }
   }

   void  Timestepper::updateTimeMatrix(SparseMatrix& oldTime, Equations::SharedIEvolutionEquation spEq, FieldComponents::Spectral::Id comp, const int idx, bool isLhs)
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

      // Update old matrix with new time dependence
      int j = 0;
      for (int k=0; k<oldTime.outerSize(); ++k)
      {
         for (SparseMatrix::InnerIterator it(oldTime,k); it; ++it)
         {
            if(static_cast<size_t>(j) < spEq->timeRow(comp,idx).first.size() && static_cast<size_t>(k) == spEq->timeRow(comp,idx).first.at(j).col())
            {
               if(static_cast<size_t>(it.row()) == spEq->timeRow(comp,idx).first.at(j).row())
               {
                  it.valueRef() += timeCoeff*spEq->timeRow(comp,idx).first.at(j).value();
                  j++;
               }
            }
         }
      }

      // Safety assert to make sure all values have been updated
      assert(static_cast<size_t>(j) == spEq->timeRow(comp,idx).first.size());
   }

   void  Timestepper::updateTimeMatrix(SparseMatrixZ& oldTime, Equations::SharedIEvolutionEquation spEq, FieldComponents::Spectral::Id comp, const int idx, bool isLhs)
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

      // Update old matrix with new time dependence
      int j = 0;
      int jj = 0;
      for (int k=0; k<oldTime.outerSize(); ++k)
      {
         for (SparseMatrixZ::InnerIterator it(oldTime,k); it; ++it)
         {
            if(static_cast<size_t>(j) < spEq->timeRow(comp,idx).first.size() && static_cast<size_t>(k) == spEq->timeRow(comp,idx).first.at(j).col())
            {
               if(static_cast<size_t>(it.row()) == spEq->timeRow(comp,idx).first.at(j).row())
               {
                  it.valueRef().real() += spEq->timeRow(comp,idx).first.at(j).value();
                  j++;
               }
            }

            if(static_cast<size_t>(jj) < spEq->timeRow(comp,idx).second.size() && static_cast<size_t>(k) == spEq->timeRow(comp,idx).second.at(jj).col())
            {
               if(static_cast<size_t>(it.row()) == spEq->timeRow(comp,idx).second.at(jj).row())
               {
                  it.valueRef().imag() += spEq->timeRow(comp,idx).second.at(jj).value();
                  jj++;
               }
            }
         }
      }

      // Safety assert to make sure all values have been updated
      assert(static_cast<size_t>(j) == spEq->timeRow(comp,idx).first.size());
      assert(static_cast<size_t>(jj) == spEq->timeRow(comp,idx).second.size());
   }
}
}
