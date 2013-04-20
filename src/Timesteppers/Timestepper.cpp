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
         // Update the time dependence in matrices
         this->updateMatrices();
         DebuggerMacro_stop("Update matrices t = ", 0);

         DebuggerMacro_start("Complex solver update", 0);
         // Update solvers from complex equation steppers
         std::vector<EquationZTimestepper>::iterator   zIt;
         for(zIt = this->mEqZStepper.begin(); zIt != this->mEqZStepper.end(); ++zIt)
         {
            zIt->updateSolver();
         }
         DebuggerMacro_stop("Complex solver update t = ", 0);

         DebuggerMacro_start("Real solver update", 0);
         // Update solvers from real equation steppers
         std::vector<EquationDTimestepper>::iterator   rIt;
         for(rIt = this->mEqDStepper.begin(); rIt != this->mEqDStepper.end(); ++rIt)
         {
            rIt->updateSolver();
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

      DebuggerMacro_start("Create timesteppers", 0);
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
      DebuggerMacro_stop("Create timesteppers t = ", 0);

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
      std::vector<EquationZTimestepper>::iterator   zIt;
      for(zIt = this->mEqZStepper.begin(); zIt != this->mEqZStepper.end(); ++zIt)
      {
         zIt->initSolver();
      }
      DebuggerMacro_stop("Complex solver init t = ", 0);

      DebuggerMacro_start("Real solver init", 0);
      // Initialise solvers from real equation steppers
      std::vector<EquationDTimestepper>::iterator   rIt;
      for(rIt = this->mEqDStepper.begin(); rIt != this->mEqDStepper.end(); ++rIt)
      {
         rIt->initSolver();
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
         if(spEq->couplingInfo(comp).solverIndex() > static_cast<int>(this->mEqZStepper.size()) - 1)
         {
            this->mEqZStepper.push_back(EquationZTimestepper(spEq->couplingInfo(comp).fieldStart()));
         }

      // Equation is part of a real system
      } else
      {
         // Add equation stepper if system index does not yet exist
         if(spEq->couplingInfo(comp).solverIndex() > static_cast<int>(this->mEqDStepper.size()) - 1)
         {
            this->mEqDStepper.push_back(EquationDTimestepper(spEq->couplingInfo(comp).fieldStart()));
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
         // Initialise the equation stepper
         if(this->mEqZStepper.at(myIdx).nSystem() == 0)
         {
            // Reserve storage for matrice and initialise vectors
            this->mEqZStepper.at(myIdx).initMatrices(ImExRK3::STEPS*nSystems);

            // Initialise field storage and information
            for(int i = 0; i < nSystems; i++)
            {
               // Create RHS and solution data storage
               this->mEqZStepper.at(myIdx).addStorage(spEq->couplingInfo(comp).systemN(i), spEq->couplingInfo(comp).rhsCols(i));

               // Build time matrix
               this->buildTimeMatrix(this->mEqZStepper.at(myIdx).rTMatrix(i), spEq, comp, i);
            }
         }

         // Build the solver matrices
         for(int i = 0; i < nSystems; i++)
         {
            // Build LHS solver matrix
            this->buildSolverMatrix(this->mEqZStepper.at(myIdx).rLHSMatrix(start+i), spEq, comp, i, true);

            // Build RHS solver matrix
            this->buildSolverMatrix(this->mEqZStepper.at(myIdx).rRHSMatrix(start+i), spEq, comp, i, false);

            // Store the start row
            startRow(i) = spEq->couplingInfo(comp).fieldIndex()*spEq->couplingInfo(comp).blockN(i);
         }

         // Store storage information
         this->mEqZStepper.at(myIdx).addInformation(myId,startRow);

      // Real matrices in linear solve
      } else
      {
         // Initialise the equation stepper
         if(this->mEqDStepper.at(myIdx).nSystem() == 0)
         {
            // Reserve storage for matrice and initialise vectors
            this->mEqDStepper.at(myIdx).initMatrices(ImExRK3::STEPS*nSystems);

            // Initialise field storage and information
            for(int i = 0; i < nSystems; i++)
            {
               // Create RHS and solution data storage
               this->mEqDStepper.at(myIdx).addStorage(spEq->couplingInfo(comp).systemN(i), spEq->couplingInfo(comp).rhsCols(i));

               // Build LHS solver matrix
               this->buildTimeMatrix(this->mEqDStepper.at(myIdx).rTMatrix(i), spEq, comp, i);
            }
         }

         // Build the solver matrices
         for(int i = 0; i < nSystems; i++)
         {
            // Build LHS solver matrix
            this->buildSolverMatrix(this->mEqDStepper.at(myIdx).rLHSMatrix(start+i), spEq, comp, i, true);

            // Build RHS solver matrix
            this->buildSolverMatrix(this->mEqDStepper.at(myIdx).rRHSMatrix(start+i), spEq, comp, i, false);

            // Store the start row
            startRow(i) = spEq->couplingInfo(comp).fieldIndex()*spEq->couplingInfo(comp).blockN(i);
         }

         // Store storage information
         this->mEqDStepper.at(myIdx).addInformation(myId,startRow);
      }
   }

   void Timestepper::updateMatrices()
   {
      // Loop over all substeps of timestepper
      for(int step = 0; step < ImExRK3::STEPS; ++step)
      {
         // Compute timestep correction coefficient for LHS matrix
         MHDFloat lhsCoeff = ImExRK3::lhsT(step)*(1.0/this->mOldDt - 1.0/this->mDt);

         // Compute timestep correction coefficient for RHS matrix
         MHDFloat rhsCoeff = ImExRK3::rhsT(step)*(1.0/this->mOldDt - 1.0/this->mDt);

         // Loop over all complex timesteppers
         std::vector<EquationZTimestepper>::iterator   zIt;
         for(zIt = this->mEqZStepper.begin(); zIt != this->mEqZStepper.end(); ++zIt)
         {
            zIt->updateTimeMatrix(lhsCoeff, rhsCoeff, step);
         }

         // Loop over all real timesteppers
         std::vector<EquationDTimestepper>::iterator   rIt;
         for(rIt = this->mEqDStepper.begin(); rIt != this->mEqDStepper.end(); ++rIt)
         {
            rIt->updateTimeMatrix(lhsCoeff, rhsCoeff, step);
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
         // Get field identity
         myId = std::make_pair((*scalEqIt)->name(), comp);

         // Get index of timesteppper
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

         // Get field identity
         myId = std::make_pair((*vectEqIt)->name(), comp);

         // Get index of timestespper
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

         // Get field identity
         myId = std::make_pair((*vectEqIt)->name(), comp);

         // Get index of timestepper
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
         // Get field identity
         myId = std::make_pair((*scalEqIt)->name(), comp);

         // Get index of timestepper
         int myIdx = (*scalEqIt)->couplingInfo(comp).solverIndex();

         // Linear solve matrices are complex
         if((*scalEqIt)->couplingInfo(comp).isComplex())
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
         // Work on first component
         comp = FieldComponents::Spectral::ONE;

         // Get field identity
         myId = std::make_pair((*vectEqIt)->name(), comp);

         // Get index of timestepper
         int myIdx = (*vectEqIt)->couplingInfo(comp).solverIndex();

         // Linear solve matrices are complex
         if((*vectEqIt)->couplingInfo(comp).isComplex())
         {
            // Get timestep input for first component
            for(int i = 0; i < this->mEqZStepper.at(myIdx).nSystem(); i++)
            {
               (*vectEqIt)->timestepInput(comp, this->mEqZStepper.at(myIdx).rRHSData(i), i, this->mEqZStepper.at(myIdx).startRow(myId,i));
            }

         // Linear solve matrices are real
         } else
         {
            // Get timestep input for first component
            for(int i = 0; i < this->mEqDStepper.at(myIdx).nSystem(); i++)
            {
               (*vectEqIt)->timestepInput(comp, this->mEqDStepper.at(myIdx).rRHSData(i), i, this->mEqDStepper.at(myIdx).startRow(myId,i));
            }
         }

         // Work on second component
         comp = FieldComponents::Spectral::TWO;

         // Get field identity
         myId = std::make_pair((*vectEqIt)->name(), comp);

         // Get index of timestepper
         myIdx = (*vectEqIt)->couplingInfo(comp).solverIndex();

         // Linear solve matrices are complex
         if((*vectEqIt)->couplingInfo(comp).isComplex())
         {
            // Get timestep input for second component
            for(int i = 0; i < this->mEqZStepper.at(myIdx).nSystem(); i++)
            {
               (*vectEqIt)->timestepInput(comp, this->mEqZStepper.at(myIdx).rRHSData(i), i, this->mEqZStepper.at(myIdx).startRow(myId,i));
            }

         // Linear solve matrices are real
         } else
         {
            // Get timestep input for second component
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
      // Storage for identity
      std::pair<PhysicalNames::Id, FieldComponents::Spectral::Id> myId;

      // Storage for component
      FieldComponents::Spectral::Id comp;

      // Loop over all scalar equations
      comp = FieldComponents::Spectral::SCALAR;
      std::vector<Equations::SharedIScalarEquation>::const_iterator scalEqIt;
      for(scalEqIt = scalEq.begin(); scalEqIt < scalEq.end(); scalEqIt++)
      {
         // Get field identity
         myId = std::make_pair((*scalEqIt)->name(), comp);

         // Get index of timestepper
         int myIdx = (*scalEqIt)->couplingInfo(comp).solverIndex();

         // Linear solve matrices are complex
         if((*scalEqIt)->couplingInfo(comp).isComplex())
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
         // Work on first component
         comp = FieldComponents::Spectral::ONE;

         // Get field identity
         myId = std::make_pair((*vectEqIt)->name(), comp);

         // Get index of timestepper
         int myIdx = (*vectEqIt)->couplingInfo(comp).solverIndex();

         // Linear solve matrices are complex
         if((*vectEqIt)->couplingInfo(comp).isComplex())
         {
            // Get timestep output for first component
            for(int i = 0; i < this->mEqZStepper.at(myIdx).nSystem(); i++)
            {
               (*vectEqIt)->timestepOutput(comp, this->mEqZStepper.at(myIdx).solution(i), i, this->mEqZStepper.at(myIdx).startRow(myId,i));
            }

         // Linear solve matrices are real
         } else
         {
            // Get timestep output for first component
            for(int i = 0; i < this->mEqDStepper.at(myIdx).nSystem(); i++)
            {
               (*vectEqIt)->timestepOutput(comp, this->mEqDStepper.at(myIdx).solution(i), i, this->mEqDStepper.at(myIdx).startRow(myId,i));
            }
         }

         // Work on second component
         comp = FieldComponents::Spectral::TWO;

         // Get field identity
         myId = std::make_pair((*vectEqIt)->name(), comp);

         // Get index of timestepper
         myIdx = (*vectEqIt)->couplingInfo(comp).solverIndex();

         // Linear solve matrices are complex
         if((*vectEqIt)->couplingInfo(comp).isComplex())
         {
            // Get timestep input for second component
            for(int i = 0; i < this->mEqZStepper.at(myIdx).nSystem(); i++)
            {
               (*vectEqIt)->timestepOutput(comp, this->mEqZStepper.at(myIdx).solution(i), i, this->mEqZStepper.at(myIdx).startRow(myId,i));
            }

         // Linear solve matrices are real
         } else
         {
            // Get timestep output for second component
            for(int i = 0; i < this->mEqDStepper.at(myIdx).nSystem(); i++)
            {
               (*vectEqIt)->timestepOutput(comp, this->mEqDStepper.at(myIdx).solution(i), i, this->mEqDStepper.at(myIdx).startRow(myId,i));
            }
         }
      }
   }

   void Timestepper::buildTimeMatrix(SparseMatrix& timeMatrix, Equations::SharedIEvolutionEquation spEq, FieldComponents::Spectral::Id comp, const int idx)
   {
      if(timeMatrix.size() == 0)
      {
         timeMatrix.resize(spEq->couplingInfo(comp).systemN(idx), spEq->couplingInfo(comp).systemN(idx));
      }

      timeMatrix += spEq->timeRow(comp, idx).first;
   }

   void Timestepper::buildTimeMatrix(SparseMatrixZ& timeMatrix, Equations::SharedIEvolutionEquation spEq, FieldComponents::Spectral::Id comp, const int idx)
   {
      if(timeMatrix.size() == 0)
      {
         timeMatrix.resize(spEq->couplingInfo(comp).systemN(idx), spEq->couplingInfo(comp).systemN(idx));
      }

      timeMatrix += spEq->timeRow(comp, idx).first.cast<MHDComplex>() + MathConstants::cI*spEq->timeRow(comp, idx).second;
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

      if(solverMatrix.size() == 0)
      {
         solverMatrix.resize(spEq->couplingInfo(comp).systemN(idx), spEq->couplingInfo(comp).systemN(idx));
      }

      // Add boundary row for LHS operator
      if(isLhs)
      {
         solverMatrix += spEq->boundaryRow(comp, idx).first;
      }

      solverMatrix += linearCoeff*spEq->linearRow(comp, idx).first - timeCoeff*spEq->timeRow(comp, idx).first;
   }

   void Timestepper::buildSolverMatrix(SparseMatrixZ& solverMatrix, Equations::SharedIEvolutionEquation spEq, FieldComponents::Spectral::Id comp, const int idx, const bool isLhs)
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

      if(solverMatrix.size() == 0)
      {
         solverMatrix.resize(spEq->couplingInfo(comp).systemN(idx), spEq->couplingInfo(comp).systemN(idx));
      }

      // Add boundary row for LHS operator
      if(isLhs)
      {
         solverMatrix += spEq->boundaryRow(comp, idx).first.cast<MHDComplex>() + MathConstants::cI*spEq->boundaryRow(comp, idx).second;
      }

      solverMatrix += linearCoeff*spEq->linearRow(comp, idx).first.cast<MHDComplex>() + MathConstants::cI*linearCoeff*spEq->linearRow(comp, idx).second - timeCoeff*spEq->timeRow(comp, idx).first.cast<MHDComplex>() - MathConstants::cI*timeCoeff*spEq->timeRow(comp, idx).second;
   }

}
}
