/** \file Timestepper.cpp
 *  \brief Implementation of a general timestepper structure
 */

// System includes
//

// External includes
//
#include <Eigen/Sparse>
#include <Eigen/KroneckerProduct>
#include <Eigen/IterativeLinearSolvers>

// Class include
//
#include "Timesteppers/Timestepper.hpp"

// Project includes
//
#include "Timesteppers/ImExRK3.hpp"

namespace GeoMHDiSCC {

   Timestepper::Timestepper()
      : mStep(0), mDt(1e-4), mTime(0.0)
   {
   }

   Timestepper::~Timestepper()
   {
   }

   bool Timestepper::finishedStep() const
   {
      return this->mStep == 0;
   }

   void Timestepper::getEqStepperType(SharedEvolutionEquation spEq, FieldComponents::Spectral::Id comp, std::map<Timestepper::FieldIdType, int>& coupled, std::vector<bool>& typeInfo)
   {
      // Get Id for current unknown field
      FieldIdType myId = std::make_pair(spEq->name(),comp);

      int pos;

      // Equation is not coupled (or is first of coupled equations)
      if(coupled.count(myId) == 0)
      {
         // Get the range of coupled fields
         CouplingInformation::FieldIteratorType it;
         CouplingInformation::FieldRangeType range = spEq->couplingInfo().fields(comp);

         pos = typeInfo.size();

         // Get current position in the equation stepper list
         for(it = range.first; it != range.second; it++)
         {
            coupled.insert(std::make_pair(it->second, pos));
         }

         typeInfo.push_back(spEq->isComplex(comp));

      // Coupled equation
      } else
      {
         // Get index of corresponding information
         pos = coupled.find(myId)->second;

         // Correct complex flag if required
         typeInfo.at(pos) = typeInfo.at(pos) || spEq->isComplex(comp);

         // Remove equation from coupled equation list
         coupled.erase(myId);
      }

      // Store preliminary equation stepper information
      this->mEqInformation.insert(std::make_pair(myId,std::make_pair(spEq->isComplex(comp), pos)));
   }

   void Timestepper::createEqStepper(SharedEvolutionEquation spEq, FieldComponents::Spectral::Id comp, std::map<Timestepper::FieldIdType, std::pair<bool,int> >& coupled, const std::vector<bool>& typeInfo)
   {
      // Get Id for current unknown field
      FieldIdType myId = std::make_pair(spEq->name(),comp);

      // Get the number of coupled fields
      int nC = spEq->couplingInfo().nFields(comp);

      // Equation is not coupled (or is first of coupled equations)
      if(coupled.count(myId) == 0)
      {
         int idx;

         // Get the range of coupled fields
         CouplingInformation::FieldIteratorType it;
         CouplingInformation::FieldRangeType range = spEq->couplingInfo().fields(comp);

         // Update data type information
         this->mEqInformation.find(myId)->second.first = typeInfo.at(this->mEqInformation.find(myId)->second.second);

         // Create correct equation timestepper
         if(this->mEqInformation.find(myId)->second.first)
         {
            idx = this->mEqZStepper.size();
            this->mEqZStepper.push_back(EquationZTimestepper(nC+1));
         } else
         {
            idx = this->mEqDStepper.size();
            this->mEqDStepper.push_back(EquationDTimestepper(nC+1));
         }

         // Update index information
         this->mEqInformation.find(myId)->second.second = idx;

         // Create list of coupled fields
         for(it = range.first; it != range.second; it++)
         {
            coupled.insert(std::make_pair(it->second, std::make_pair(this->mEqInformation.find(myId)->second.first,idx)));
         }

      // Coupled equation
      } else
      {
         // Update data type information
         this->mEqInformation.find(myId)->second.first = coupled.find(myId)->second.first;

         // Update index information
         this->mEqInformation.find(myId)->second.second = coupled.find(myId)->second.second;

         // Remove equation from coupled equation list
         coupled.erase(myId);
      }
   }

   void Timestepper::init(const std::vector<SharedScalarEquation>& scalEq, const std::vector<SharedVectorEquation>& vectEq)
   {
      // Storage for the field coupling information
      std::map<FieldIdType, int>  coupled;
      std::vector<bool>  typeInfo;

      // Loop over all scalar equations
      std::vector<SharedScalarEquation>::const_iterator scalEqIt;
      for(scalEqIt = scalEq.begin(); scalEqIt < scalEq.end(); scalEqIt++)
      {
         // Get type information for the equation steppers
         this->getEqStepperType((*scalEqIt), FieldComponents::Spectral::NONE, coupled, typeInfo);
      }

      // Loop over all vector equations
      std::vector<SharedVectorEquation>::const_iterator vectEqIt;
      for(vectEqIt = vectEq.begin(); vectEqIt < vectEq.end(); vectEqIt++)
      {
         // Get type information for the equation steppers for the toroidal component
         this->getEqStepperType((*vectEqIt), FieldComponents::Spectral::TOROIDAL, coupled, typeInfo);

         // Get type information for the equation steppers for the poloidal component
         this->getEqStepperType((*vectEqIt), FieldComponents::Spectral::POLOIDAL, coupled, typeInfo);
      }

      // Storage for the field coupling information
      std::map<FieldIdType, std::pair<bool,int> >  coupledType;

      // Loop over all scalar equations
      for(scalEqIt = scalEq.begin(); scalEqIt < scalEq.end(); scalEqIt++)
      {
         // Get type information for the equation steppers
         this->createEqStepper((*scalEqIt), FieldComponents::Spectral::NONE, coupledType, typeInfo);
      }

      // Loop over all vector equations
      for(vectEqIt = vectEq.begin(); vectEqIt < vectEq.end(); vectEqIt++)
      {
         // Get type information for the equation steppers for the toroidal component
         this->createEqStepper((*vectEqIt), FieldComponents::Spectral::TOROIDAL, coupledType, typeInfo);

         // Get type information for the equation steppers for the poloidal component
         this->createEqStepper((*vectEqIt), FieldComponents::Spectral::POLOIDAL, coupledType, typeInfo);
      }

      // Loop over all substeps of timestepper
      for(this->mStep = 0; this->mStep < ImExRK3::STEPS; this->mStep++)
      {
         // Loop over all scalar equations
         for(scalEqIt = scalEq.begin(); scalEqIt < scalEq.end(); scalEqIt++)
         {
            // Create (coupled) matrices
            this->createMatrices((*scalEqIt), FieldComponents::Spectral::NONE);
         }

         // Loop over all vector equations
         for(vectEqIt = vectEq.begin(); vectEqIt < vectEq.end(); vectEqIt++)
         {
            // Create (coupled) matrices
            this->createMatrices((*vectEqIt), FieldComponents::Spectral::TOROIDAL);

            // Create (coupled) matrices
            this->createMatrices((*vectEqIt), FieldComponents::Spectral::POLOIDAL);
         }
      }

      // Initialise solvers from equation steppers
      for(int i = 0; i < this->mEqZStepper.size(); i++)
      {
         this->mEqZStepper.at(i).initSolver();
      }

      // Initialise solvers from real equation steppers
      for(int i = 0; i < this->mEqDStepper.size(); i++)
      {
         this->mEqDStepper.at(i).initSolver();
      }

      // Reset the step index
      this->mStep = 0;

      // Initialise with initial state
      this->initSolution(scalEq, vectEq);
   }

   void Timestepper::createMatrices(SharedEvolutionEquation spEq, FieldComponents::Spectral::Id comp)
   {
      // ID of the current field
      FieldIdType myId = std::make_pair(spEq->name(),comp);

      // Get internal information
      std::pair<int,ArrayI>  interInfo = spEq->couplingInfo().internal(comp);

      // Get equation stepper information
      std::pair<bool,int>  eqInfo = this->mEqInformation.find(myId)->second;

      // Get the number of coupled fields
      int nC = spEq->couplingInfo().nFields(comp);

      // start index for matrices
      int start = this->mStep*interInfo.first;

      // Start row for storage information
      ArrayI startRow(interInfo.first);

      // Storage for the field index
      int fieldIdx;

      // Complex matrices in linear solve
      if(eqInfo.first)
      {
         fieldIdx = this->mEqZStepper.at(eqInfo.second).current();

         for(int i = 0; i < interInfo.first; i++)
         {
            if(fieldIdx == 0)
            {
               // Complete RHS matrix
               this->mEqZStepper.at(eqInfo.second).addLHSMatrix(this->buildLHSMatrix(spEq,comp,i,nC,fieldIdx));

               // Complete RHS matrix
               this->mEqZStepper.at(eqInfo.second).addRHSMatrix(this->buildRHSMatrix(spEq,comp,i,nC,fieldIdx));

               if(this->mStep == 0)
               {
                  // Create RHS and solution data storage
                  this->mEqZStepper.at(eqInfo.second).addStorage(interInfo.second(i));
               }

               // Set the start row
               startRow(i) = 0;
            } else
            {
               // Complete RHS matrix
               this->mEqZStepper.at(eqInfo.second).completeLHSMatrix(start+i,this->buildLHSMatrix(spEq,comp,i,nC,fieldIdx));

               // Complete RHS matrix
               this->mEqZStepper.at(eqInfo.second).completeRHSMatrix(start+i,this->buildRHSMatrix(spEq,comp,i,nC,fieldIdx));

               // Set the start row
               startRow(i) = spEq->rowShift(comp, i);
            }
         }

         if(this->mStep == 0)
         {
            // Store storage information
            this->mEqZStepper.at(eqInfo.second).addInformation(myId,startRow);
         }

         // Move field counter to next field
         this->mEqZStepper.at(eqInfo.second).next();

      // Real matrices in linear solve
      } else
      {
         fieldIdx = this->mEqDStepper.at(eqInfo.second).current();

         for(int i = 0; i < interInfo.first; i++)
         {
            if(fieldIdx == 0)
            {
               // Complete RHS matrix
               this->mEqDStepper.at(eqInfo.second).addLHSMatrix(this->buildLHSMatrix(spEq,comp,i,nC,fieldIdx));

               // Complete RHS matrix
               this->mEqDStepper.at(eqInfo.second).addRHSMatrix(this->buildRHSMatrix(spEq,comp,i,nC,fieldIdx));

               if(this->mStep == 0)
               {
                  // Create RHS and solution data storage
                  this->mEqDStepper.at(eqInfo.second).addStorage(interInfo.second(i));
               }

               // Set the start row
               startRow(i) = 0;
            } else
            {
               // Complete RHS matrix
               this->mEqDStepper.at(eqInfo.second).completeLHSMatrix(start+i,this->buildLHSMatrix(spEq,comp,i,nC,fieldIdx));

               // Complete RHS matrix
               this->mEqDStepper.at(eqInfo.second).completeRHSMatrix(start+i,this->buildRHSMatrix(spEq,comp,i,nC,fieldIdx));

               // Set the start row
               startRow(i) = spEq->rowShift(comp, i);
            }
         }

         if(this->mStep == 0)
         {
            // Store storage information
            this->mEqDStepper.at(eqInfo.second).addInformation(myId,startRow);
         }

         // Move field counter to next field
         this->mEqDStepper.at(eqInfo.second).next();
      }
   }

   void Timestepper::update()
   {
      this->mTime += this->mDt;
      std::cerr << this->mTime << std::endl;
   }

   void Timestepper::initSolution(const std::vector<SharedScalarEquation>& scalEq, const std::vector<SharedVectorEquation>& vectEq)
   {
      // Storage for information and identity
      std::pair<bool,int> eqInfo;
      FieldIdType myId;

      // Loop over all scalar equations
      std::vector<SharedScalarEquation>::const_iterator scalEqIt;
      for(scalEqIt = scalEq.begin(); scalEqIt < scalEq.end(); scalEqIt++)
      {
         // Get identity and corresponding equation information
         myId = std::make_pair((*scalEqIt)->name(), FieldComponents::Spectral::NONE);
         // Get equation stepper information
         eqInfo = this->mEqInformation.find(myId)->second;

         // Linear solve matrices are complex
         if(eqInfo.first)
         {
            // Get timestep input
            for(int i = 0; i < this->mEqZStepper.at(eqInfo.second).nSystem(); i++)
            {
               (*scalEqIt)->copyTInput(FieldComponents::Spectral::NONE, this->mEqZStepper.at(eqInfo.second).rSolution(i), i, this->mEqZStepper.at(eqInfo.second).startRow(myId,i));
            }

         // Linear solve matrices are real
         } else
         {
            // Get timestep input
            for(int i = 0; i < this->mEqDStepper.at(eqInfo.second).nSystem(); i++)
            {
               (*scalEqIt)->copyTInput(FieldComponents::Spectral::NONE, this->mEqDStepper.at(eqInfo.second).rSolution(i), i, this->mEqDStepper.at(eqInfo.second).startRow(myId,i));
            }
         }
      }

      // Loop over all vector equations
      std::vector<SharedVectorEquation>::const_iterator vectEqIt;
      for(vectEqIt = vectEq.begin(); vectEqIt < vectEq.end(); vectEqIt++)
      {
         // Get identity and corresponding equation information for toroidal component
         myId = std::make_pair((*vectEqIt)->name(), FieldComponents::Spectral::TOROIDAL);
         // Get equation stepper information
         eqInfo = this->mEqInformation.find(myId)->second;

         // Linear solve matrices are complex
         if(eqInfo.first)
         {
            // Get timestep input for toroidal component
            for(int i = 0; i < this->mEqZStepper.at(eqInfo.second).nSystem(); i++)
            {
               (*vectEqIt)->copyTInput(FieldComponents::Spectral::TOROIDAL, this->mEqZStepper.at(eqInfo.second).rSolution(i), i, this->mEqZStepper.at(eqInfo.second).startRow(myId,i));
            }

         // Linear solve matrices are real
         } else
         {
            // Get timestep input for toroidal component
            for(int i = 0; i < this->mEqDStepper.at(eqInfo.second).nSystem(); i++)
            {
               (*vectEqIt)->copyTInput(FieldComponents::Spectral::TOROIDAL, this->mEqDStepper.at(eqInfo.second).rSolution(i), i, this->mEqDStepper.at(eqInfo.second).startRow(myId,i));
            }
         }

         // Get identity and corresponding equation information for poloidal component
         myId = std::make_pair((*vectEqIt)->name(), FieldComponents::Spectral::POLOIDAL);
         // Get equation stepper information
         eqInfo = this->mEqInformation.find(myId)->second;

         // Linear solve matrices are complex
         if(eqInfo.first)
         {
            // Get timestep input for poloidal component
            for(int i = 0; i < this->mEqZStepper.at(eqInfo.second).nSystem(); i++)
            {
               (*vectEqIt)->copyTInput(FieldComponents::Spectral::POLOIDAL, this->mEqZStepper.at(eqInfo.second).rSolution(i), i, this->mEqZStepper.at(eqInfo.second).startRow(myId,i));
            }

         // Linear solve matrices are real
         } else
         {
            // Get timestep input for poloidal component
            for(int i = 0; i < this->mEqDStepper.at(eqInfo.second).nSystem(); i++)
            {
               (*vectEqIt)->copyTInput(FieldComponents::Spectral::POLOIDAL, this->mEqDStepper.at(eqInfo.second).rSolution(i), i, this->mEqDStepper.at(eqInfo.second).startRow(myId,i));
            }
         }
      }
   }

   void Timestepper::stepForward(const std::vector<SharedScalarEquation>& scalEq, const std::vector<SharedVectorEquation>& vectEq)
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

   void Timestepper::getInput(const std::vector<SharedScalarEquation>& scalEq, const std::vector<SharedVectorEquation>& vectEq)
   {
      // Storage for information and identity
      std::pair<bool,int> eqInfo;
      FieldIdType myId;

      // Loop over all scalar equations
      std::vector<SharedScalarEquation>::const_iterator scalEqIt;
      for(scalEqIt = scalEq.begin(); scalEqIt < scalEq.end(); scalEqIt++)
      {
         // Get identity and corresponding equation information
         myId = std::make_pair((*scalEqIt)->name(), FieldComponents::Spectral::NONE);
         // Get equation stepper information
         eqInfo = this->mEqInformation.find(myId)->second;

         // Linear solve matrices are complex
         if(eqInfo.first)
         {
            // Get timestep input
            for(int i = 0; i < this->mEqZStepper.at(eqInfo.second).nSystem(); i++)
            {
               (*scalEqIt)->timestepInput(FieldComponents::Spectral::NONE, this->mEqZStepper.at(eqInfo.second).rRHSData(i), i, this->mEqZStepper.at(eqInfo.second).startRow(myId,i));
            }

         // Linear solve matrices are real
         } else
         {
            // Get timestep input
            for(int i = 0; i < this->mEqDStepper.at(eqInfo.second).nSystem(); i++)
            {
               (*scalEqIt)->timestepInput(FieldComponents::Spectral::NONE, this->mEqDStepper.at(eqInfo.second).rRHSData(i), i, this->mEqDStepper.at(eqInfo.second).startRow(myId,i));
            }
         }
      }

      // Loop over all vector equations
      std::vector<SharedVectorEquation>::const_iterator vectEqIt;
      for(vectEqIt = vectEq.begin(); vectEqIt < vectEq.end(); vectEqIt++)
      {
         // Get identity and corresponding equation information for toroidal component
         myId = std::make_pair((*vectEqIt)->name(), FieldComponents::Spectral::TOROIDAL);
         // Get equation stepper information
         eqInfo = this->mEqInformation.find(myId)->second;

         // Linear solve matrices are complex
         if(eqInfo.first)
         {
            // Get timestep input for toroidal component
            for(int i = 0; i < this->mEqZStepper.at(eqInfo.second).nSystem(); i++)
            {
               (*vectEqIt)->timestepInput(FieldComponents::Spectral::TOROIDAL, this->mEqZStepper.at(eqInfo.second).rRHSData(i), i, this->mEqZStepper.at(eqInfo.second).startRow(myId,i));
            }

         // Linear solve matrices are real
         } else
         {
            // Get timestep input for toroidal component
            for(int i = 0; i < this->mEqDStepper.at(eqInfo.second).nSystem(); i++)
            {
               (*vectEqIt)->timestepInput(FieldComponents::Spectral::TOROIDAL, this->mEqDStepper.at(eqInfo.second).rRHSData(i), i, this->mEqDStepper.at(eqInfo.second).startRow(myId,i));
            }
         }

         // Get identity and corresponding equation information for poloidal component
         myId = std::make_pair((*vectEqIt)->name(), FieldComponents::Spectral::POLOIDAL);
         // Get equation stepper information
         eqInfo = this->mEqInformation.find(myId)->second;

         // Linear solve matrices are complex
         if(eqInfo.first)
         {
            // Get timestep input for poloidal component
            for(int i = 0; i < this->mEqZStepper.at(eqInfo.second).nSystem(); i++)
            {
               (*vectEqIt)->timestepInput(FieldComponents::Spectral::POLOIDAL, this->mEqZStepper.at(eqInfo.second).rRHSData(i), i, this->mEqZStepper.at(eqInfo.second).startRow(myId,i));
            }

         // Linear solve matrices are real
         } else
         {
            // Get timestep input for poloidal component
            for(int i = 0; i < this->mEqDStepper.at(eqInfo.second).nSystem(); i++)
            {
               (*vectEqIt)->timestepInput(FieldComponents::Spectral::POLOIDAL, this->mEqDStepper.at(eqInfo.second).rRHSData(i), i, this->mEqDStepper.at(eqInfo.second).startRow(myId,i));
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
      // Compute RHS component for complex linear systems
      std::vector<EquationZTimestepper>::iterator   zIt;
      for(zIt = this->mEqZStepper.begin(); zIt != this->mEqZStepper.end(); ++zIt)
      {
         // Compute linear solve RHS
         zIt->solve(this->mStep);
      }

      // Compute RHS component for real linear systems
      std::vector<EquationDTimestepper>::iterator   dIt;
      for(dIt = this->mEqDStepper.begin(); dIt != this->mEqDStepper.end(); ++dIt)
      {
         // Compute linear solve RHS
         dIt->solve(this->mStep);
      }
   }

   void Timestepper::transferOutput(const std::vector<SharedScalarEquation>& scalEq, const std::vector<SharedVectorEquation>& vectEq)
   {
      // Storage for information and identity
      std::pair<bool,int> eqInfo;
      FieldIdType myId;

      // Loop over all scalar equations
      std::vector<SharedScalarEquation>::const_iterator scalEqIt;
      for(scalEqIt = scalEq.begin(); scalEqIt < scalEq.end(); scalEqIt++)
      {
         // Get identity and corresponding equation information
         myId = std::make_pair((*scalEqIt)->name(), FieldComponents::Spectral::NONE);
         // Get equation stepper information
         eqInfo = this->mEqInformation.find(myId)->second;

         // Linear solve matrices are complex
         if(eqInfo.first)
         {
            // Get timestep output
            for(int i = 0; i < this->mEqZStepper.at(eqInfo.second).nSystem(); i++)
            {
               (*scalEqIt)->timestepOutput(FieldComponents::Spectral::NONE, this->mEqZStepper.at(eqInfo.second).solution(i), i, this->mEqZStepper.at(eqInfo.second).startRow(myId,i));
            }

         // Linear solve matrices are real
         } else
         {
            // Get timestep output
            for(int i = 0; i < this->mEqDStepper.at(eqInfo.second).nSystem(); i++)
            {
               (*scalEqIt)->timestepOutput(FieldComponents::Spectral::NONE, this->mEqDStepper.at(eqInfo.second).solution(i), i, this->mEqDStepper.at(eqInfo.second).startRow(myId,i));
            }
         }
      }

      // Loop over all vector equations
      std::vector<SharedVectorEquation>::const_iterator vectEqIt;
      for(vectEqIt = vectEq.begin(); vectEqIt < vectEq.end(); vectEqIt++)
      {
         // Get identity and corresponding equation information for toroidal component
         myId = std::make_pair((*vectEqIt)->name(), FieldComponents::Spectral::TOROIDAL);
         // Get equation stepper information
         eqInfo = this->mEqInformation.find(myId)->second;

         // Linear solve matrices are complex
         if(eqInfo.first)
         {
            // Get timestep output for toroidal component
            for(int i = 0; i < this->mEqZStepper.at(eqInfo.second).nSystem(); i++)
            {
               (*vectEqIt)->timestepOutput(FieldComponents::Spectral::TOROIDAL, this->mEqZStepper.at(eqInfo.second).solution(i), i, this->mEqZStepper.at(eqInfo.second).startRow(myId,i));
            }

         // Linear solve matrices are real
         } else
         {
            // Get timestep output for toroidal component
            for(int i = 0; i < this->mEqDStepper.at(eqInfo.second).nSystem(); i++)
            {
               (*vectEqIt)->timestepOutput(FieldComponents::Spectral::TOROIDAL, this->mEqDStepper.at(eqInfo.second).solution(i), i, this->mEqDStepper.at(eqInfo.second).startRow(myId,i));
            }
         }

         // Get identity and corresponding equation information for poloidal component
         myId = std::make_pair((*vectEqIt)->name(), FieldComponents::Spectral::POLOIDAL);
         // Get equation stepper information
         eqInfo = this->mEqInformation.find(myId)->second;

         // Linear solve matrices are complex
         if(eqInfo.first)
         {
            // Get timestep input for poloidal component
            for(int i = 0; i < this->mEqZStepper.at(eqInfo.second).nSystem(); i++)
            {
               (*vectEqIt)->timestepOutput(FieldComponents::Spectral::POLOIDAL, this->mEqZStepper.at(eqInfo.second).solution(i), i, this->mEqZStepper.at(eqInfo.second).startRow(myId,i));
            }

         // Linear solve matrices are real
         } else
         {
            // Get timestep output for poloidal component
            for(int i = 0; i < this->mEqDStepper.at(eqInfo.second).nSystem(); i++)
            {
               (*vectEqIt)->timestepOutput(FieldComponents::Spectral::POLOIDAL, this->mEqDStepper.at(eqInfo.second).solution(i), i, this->mEqDStepper.at(eqInfo.second).startRow(myId,i));
            }
         }
      }
   }

   DecoupledZSparse  Timestepper::buildLHSMatrix(SharedEvolutionEquation spEq, FieldComponents::Spectral::Id comp, const int idx, const int nC, const int cRow)
   {
      // Matrix coefficients
      MHDFloat timeCoeff = ImExRK3::lhsT(this->mStep)*1.0/this->mDt;
      MHDFloat linearCoeff = ImExRK3::lhsL(this->mStep);

      DecoupledZSparse   lhs;

      // Addition storage for coupled case
      SparseMatrix   coupledT;
      SparseMatrix   coupledC;

      // Create real LHS matrix
      lhs.first =  linearCoeff*spEq->linearMatrix(comp,idx).first - timeCoeff*spEq->timeMatrix(comp,idx).first;
      // Create imaginary LHS matrix
      lhs.second = linearCoeff*spEq->linearMatrix(comp,idx).second - timeCoeff*spEq->timeMatrix(comp,idx).second;

      // Create coupled matrix
      if(nC > 0)
      {
         // Block selection matrix for coupled system
         SparseMatrix  blockT(nC+1,nC+1);
         SparseMatrix  blockC(nC+1,nC+1);

         // Single element is required
         blockT.reserve(1);
         blockC.reserve(1);

         // Create matrix blocks
         blockT.insert(cRow,cRow) = 1.0;
         blockC.insert(cRow,1-cRow) = linearCoeff;
         
         // Create real coupled LHS matrix
         Eigen::kroneckerProduct(blockT,lhs.first, coupledT);
         Eigen::kroneckerProduct(blockC,spEq->couplingMatrix(comp,idx).first, coupledC);
         lhs.first = coupledT + coupledC;

         // Create imaginary coupled LHS matrix
         Eigen::kroneckerProduct(blockT,lhs.second, coupledT);
         Eigen::kroneckerProduct(blockC,spEq->couplingMatrix(comp,idx).second, coupledC);
         lhs.second = coupledT + coupledC;

         // Add boundary conditions
         //
         // Make unit blocks (don't want rescaled boundary conditions)
         blockT.coeffRef(cRow,cRow) = 1.0;
         blockC.coeffRef(cRow,1-cRow) = 1.0;

         // Real coupled boundary conditions
         Eigen::kroneckerProduct(blockT,spEq->bcMatrix(comp,idx).first, coupledT);
         Eigen::kroneckerProduct(blockC,spEq->cbcMatrix(comp,idx).first, coupledC);
         lhs.first += coupledT + coupledC;

         // Imaginary coupled boundary conditions
         Eigen::kroneckerProduct(blockT,spEq->bcMatrix(comp,idx).second, coupledT);
         Eigen::kroneckerProduct(blockC,spEq->cbcMatrix(comp,idx).second, coupledC);
         lhs.second += coupledT + coupledC;

      // Add boundary conditions for uncoupled system
      } else
      {
         // Real boundary conditions
         lhs.first += spEq->bcMatrix(comp,idx).first;

         // Imaginary boundary conditions
         lhs.second += spEq->bcMatrix(comp,idx).second;
      }

      return lhs;
   }

   DecoupledZSparse  Timestepper::buildRHSMatrix(SharedEvolutionEquation spEq, FieldComponents::Spectral::Id comp, const int idx, const int nC, const int cRow)
   {
      // Matrix coefficients
      MHDFloat timeCoeff = -ImExRK3::rhsT(this->mStep)*1.0/this->mDt;
      MHDFloat linearCoeff = -ImExRK3::rhsL(this->mStep);

      DecoupledZSparse   rhs;

      // Addition storage for coupled case
      SparseMatrix   coupledT;
      SparseMatrix   coupledC;

      // Create real LHS matrix
      rhs.first = linearCoeff*spEq->linearMatrix(comp,idx).first + timeCoeff*spEq->timeMatrix(comp,idx).first;
      // Create imaginary LHS matrix
      rhs.second = linearCoeff*spEq->linearMatrix(comp,idx).second + timeCoeff*spEq->timeMatrix(comp,idx).second;

      // Create coupled matrix
      if(nC > 0)
      {
         // Block selection matrix for coupled system
         SparseMatrix  blockT(nC+1,nC+1);
         SparseMatrix  blockC(nC+1,nC+1);

         // Single element is required
         blockT.reserve(1);
         blockC.reserve(1);

         // Create correct equation timestepper
         blockT.insert(cRow,cRow) = 1.0;
         blockC.insert(cRow,1-cRow) = linearCoeff;

         Eigen::kroneckerProduct(blockT,rhs.first, coupledT);
         Eigen::kroneckerProduct(blockC,spEq->couplingMatrix(comp,idx).first, coupledC);
         rhs.first = coupledT + coupledC;

         Eigen::kroneckerProduct(blockT,rhs.second, coupledT);
         Eigen::kroneckerProduct(blockC,spEq->couplingMatrix(comp,idx).second, coupledC);
         rhs.second = coupledT + coupledC;
      }

      return rhs;
   }
}
