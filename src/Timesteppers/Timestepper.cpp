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

namespace Timestep {

   Timestepper::Timestepper()
      : mStep(0), mDt(1e-3), mTime(0.0)
   {
   }

   Timestepper::~Timestepper()
   {
   }

   bool Timestepper::finishedStep() const
   {
      return this->mStep == 0;
   }

   void Timestepper::update()
   {
      this->mTime += this->mDt;
   }

   void Timestepper::changeTimestep(MHDFloat dt)
   {
      // Safety assert
      assert(dt > 0);

      // Set the new timestep
      this->mDt = dt;

      // Update the timestep matrices
      ///  \mhdBug Does not yet update the matrices!!!
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

   void Timestepper::init(const std::vector<Equations::SharedIScalarEquation>& scalEq, const std::vector<Equations::SharedIVectorEquation>& vectEq)
   {
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

      //
      // Initialise the solvers and the initial state
      //

      // Initialise solvers from complex equation steppers
      for(size_t i = 0; i < this->mEqZStepper.size(); i++)
      {
         this->mEqZStepper.at(i).initSolver();
      }

      // Initialise solvers from real equation steppers
      for(size_t i = 0; i < this->mEqDStepper.size(); i++)
      {
         this->mEqDStepper.at(i).initSolver();
      }

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

      // Equation is not coupled or is the first of coupled set
      } else
      {
         // Get a new index for the system
         int idx = this->mTimeCoupling.newIndex();

         // Add equation unknown to time coupling information
         this->mTimeCoupling.addField(myId, spEq->isComplex(comp), idx);

         // Get the range of coupled fields
         Equations::CouplingInformation::field_iterator it;
         Equations::CouplingInformation::field_iterator_range range = spEq->couplingInfo().fieldRange(comp);

         // Loop over coupled fields and add them to timestep coupling information
         for(it = range.first; it != range.second; ++it)
         {
            this->mTimeCoupling.addField(it->second, false, idx);
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
            this->mEqZStepper.push_back(EquationZTimestepper(nC+1));
         } else
         {
            idx = this->mEqDStepper.size();
            this->mEqDStepper.push_back(EquationDTimestepper(nC+1));
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

         for(int i = 0; i < interInfo.first; i++)
         {
            if(fieldIdx == 0)
            {
               // Complete RHS matrix
               this->mEqZStepper.at(myIdx).addLHSMatrix(this->buildLHSMatrix(spEq,comp,i,nC,fieldIdx));

               // Complete RHS matrix
               this->mEqZStepper.at(myIdx).addRHSMatrix(this->buildRHSMatrix(spEq,comp,i,nC,fieldIdx));

               if(this->mStep == 0)
               {
                  // Create RHS and solution data storage
                  this->mEqZStepper.at(myIdx).addStorage(interInfo.second(i));
               }

               // Set the start row
               startRow(i) = 0;
            } else
            {
               // Complete RHS matrix
               this->mEqZStepper.at(myIdx).completeLHSMatrix(start+i,this->buildLHSMatrix(spEq,comp,i,nC,fieldIdx));

               // Complete RHS matrix
               this->mEqZStepper.at(myIdx).completeRHSMatrix(start+i,this->buildRHSMatrix(spEq,comp,i,nC,fieldIdx));

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

         for(int i = 0; i < interInfo.first; i++)
         {
            if(fieldIdx == 0)
            {
               // Complete RHS matrix
               this->mEqDStepper.at(myIdx).addLHSMatrix(this->buildLHSMatrix(spEq,comp,i,nC,fieldIdx));

               // Complete RHS matrix
               this->mEqDStepper.at(myIdx).addRHSMatrix(this->buildRHSMatrix(spEq,comp,i,nC,fieldIdx));

               if(this->mStep == 0)
               {
                  // Create RHS and solution data storage
                  this->mEqDStepper.at(myIdx).addStorage(interInfo.second(i));
               }

               // Set the start row
               startRow(i) = 0;
            } else
            {
               // Complete RHS matrix
               this->mEqDStepper.at(myIdx).completeLHSMatrix(start+i,this->buildLHSMatrix(spEq,comp,i,nC,fieldIdx));

               // Complete RHS matrix
               this->mEqDStepper.at(myIdx).completeRHSMatrix(start+i,this->buildRHSMatrix(spEq,comp,i,nC,fieldIdx));

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

   DecoupledZSparse  Timestepper::buildLHSMatrix(Equations::SharedIEvolutionEquation spEq, FieldComponents::Spectral::Id comp, const int idx, const int nC, const int cRow)
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

   DecoupledZSparse  Timestepper::buildRHSMatrix(Equations::SharedIEvolutionEquation spEq, FieldComponents::Spectral::Id comp, const int idx, const int nC, const int cRow)
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
}
