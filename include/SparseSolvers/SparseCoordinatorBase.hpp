/** 
 * @file SparseCoordinatorBase.hpp
 * @brief Implementation of the base for a general sparse solver coordinator
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef SPARSECOORDINATORBASE_HPP
#define SPARSECOORDINATORBASE_HPP

// Configuration includes
//
#include "SmartPointers/SharedPtrMacro.h"

// System includes
//

// External includes
//

// Project includes
//
#include "SparseSolvers/SparseLinearSolver.hpp"
#include "Equations/IScalarEquation.hpp"
#include "Equations/IVectorEquation.hpp"

namespace GeoMHDiSCC {

namespace Solver {

   /**
    * @brief Implementation of the base for a general sparse solver coordinator
    */
   template <typename TSharedZSolver, typename TSharedRZSolver> class SparseCoordinatorBase
   {
      public:
         /// Typedef for an iterator to a complex linear solver
         typedef typename std::vector<TSharedZSolver>::iterator   SolverZ_iterator;

         /// Typedef for an iterator to a real linear solver
         typedef typename std::vector<TSharedRZSolver>::iterator   SolverRZ_iterator;

         /// Typedef for a shared scalar equation iterator
         typedef std::vector<Equations::SharedIScalarEquation>::iterator   ScalarEquation_iterator;

         /// Typedef for a shared vector equation iterator
         typedef std::vector<Equations::SharedIVectorEquation>::iterator   VectorEquation_iterator;

         /// Typedef for a shared scalar equation range
         typedef std::pair<ScalarEquation_iterator, ScalarEquation_iterator>  ScalarEquation_range;

         /// Typedef for a shared vector equation range
         typedef std::pair<VectorEquation_iterator, VectorEquation_iterator>  VectorEquation_range;

         /// Typedef for a shared scalar variable map
         typedef std::map<PhysicalNames::Id, Datatypes::SharedScalarVariableType>  ScalarVariable_map;

         /// Typedef for a shared vector variable map
         typedef std::map<PhysicalNames::Id, Datatypes::SharedVectorVariableType>  VectorVariable_map;

         /// Typedef for a shared scalar variable iterator
         typedef ScalarVariable_map::iterator  ScalarVariable_iterator;

         /// Typedef for a shared vector variable iterator
         typedef VectorVariable_map::iterator  VectorVariable_iterator;

         /// Typedef for a shared scalar variable range
         typedef std::pair<ScalarVariable_iterator, ScalarVariable_iterator>  ScalarVariable_range;

         /// Typedef for a shared vector variable range
         typedef std::pair<VectorVariable_iterator, VectorVariable_iterator>  VectorVariable_range;

         /**
          * @brief Constructor
          */
         SparseCoordinatorBase();

         /**
          * @brief Destructor
          */
         virtual ~SparseCoordinatorBase();

         /**
          * @brief Finished computation of solver step?
          */
         bool finishedStep() const;

         /**
          * @brief Initialise solver coordinator
          *
          * @param scalEq  Shared scalar equations
          * @param vectEq  Shared vector equations
          */
         virtual void init(const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq) = 0;
         
      protected:
         /**
          * @brief Create a real linear solver
          */
         virtual void addSolverRZ(const int start) = 0;

         /**
          * @brief Create a complex linear solver
          */
         virtual void addSolverZ(const int start) = 0;

         /**
          * @brief Get the solver input independently of solver type
          */
         template <typename TEquationIt, typename TSolverIt> void getSolverInput(const TEquationIt eqIt, const SpectralFieldId id, const TSolverIt solveIt, const ScalarVariable_map& scalVar, const VectorVariable_map& vectVar);

         /**
          * @brief Create the correct solver
          */
         void createSolver(Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp);

         /**
          * @brief Create storage for the solvers
          */
         void createStorage(Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp);

         /**
          * @brief Setup the storage and related information independently of solver type
          */
         template <typename TSolverIt> void setupStorage(Equations::SharedIEquation spEq, const SpectralFieldId id, const TSolverIt solveIt);

         /**
          * @brief Initialise the solution
          *
          * @param scalEq Scalar equations
          * @param vectEq Vector equations
          */
         void initSolution(const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq);

         /**
          * @brief Update equation input to solver
          *
          * @param scalEq Scalar equations
          * @param vectEq Vector equations
          */
         void getInput(const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq, const ScalarVariable_map& scalVar, const VectorVariable_map& vectVar);

         /**
          * @brief Update equation unkowns with solver output 
          *
          * @param scalEq Shared scalar equations
          * @param vectEq Shared vector equations
          */
         void transferOutput(const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq);

         /**
          * @brief Number of solver substeps
          */
         int   mNStep;

         /**
          * @brief Current solver step
          */
         int   mStep;

         /**
          * @brief Vector of (coupled) real solvers
          */
         std::vector<TSharedRZSolver> mRZSolvers;

         /**
          * @brief Vector of (coupled) complex solvers
          */
         std::vector<TSharedZSolver> mZSolvers;

      private:
   };

   template <typename TSharedZSolver, typename TSharedRZSolver> SparseCoordinatorBase<TSharedZSolver,TSharedRZSolver>::SparseCoordinatorBase()
      : mNStep(1), mStep(0)
   {
   }

   template <typename TSharedZSolver, typename TSharedRZSolver> SparseCoordinatorBase<TSharedZSolver,TSharedRZSolver>::~SparseCoordinatorBase()
   {
   }

   template <typename TSharedZSolver, typename TSharedRZSolver> bool SparseCoordinatorBase<TSharedZSolver,TSharedRZSolver>::finishedStep() const
   {
      return this->mStep == 0;
   }

   template <typename TSharedZSolver, typename TSharedRZSolver> void SparseCoordinatorBase<TSharedZSolver,TSharedRZSolver>::createSolver(Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp)
   {
      // Equation is part of a complex system
      if(spEq->couplingInfo(comp).isComplex())
      {
         // Add linear solver if system index does not yet exist
         if(spEq->couplingInfo(comp).solverIndex() > static_cast<int>(this->mZSolvers.size()) - 1)
         {
            this->addSolverZ(spEq->couplingInfo(comp).fieldStart());
         }

      // Equation is part of a real system
      } else
      {
         // Add equation stepper if system index does not yet exist
         if(spEq->couplingInfo(comp).solverIndex() > static_cast<int>(this->mRZSolvers.size()) - 1)
         {
            this->addSolverRZ(spEq->couplingInfo(comp).fieldStart());
         }
      }
   }

   template <typename TSharedZSolver, typename TSharedRZSolver> void SparseCoordinatorBase<TSharedZSolver,TSharedRZSolver>::createStorage(Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp)
   {
      // ID of the current field
      SpectralFieldId myId = std::make_pair(spEq->name(),comp);

      // Index of the current field
      int myIdx = spEq->couplingInfo(myId.second).solverIndex();

      // Complex matrices in linear solve
      if(spEq->couplingInfo(myId.second).isComplex())
      {
         // Create iterator to current complex solver
         SolverZ_iterator solZIt = this->mZSolvers.begin();
         std::advance(solZIt, myIdx);

         // setup storage and information
         this->setupStorage(spEq, myId, solZIt);

      // Real matrices in linear solve
      } else
      {
         // Create iterator to current real solver
         SolverRZ_iterator solRZIt = this->mRZSolvers.begin();
         std::advance(solRZIt, myIdx);

         // setup storage and information
         this->setupStorage(spEq, myId, solRZIt);
      }
   }

   template <typename TSharedZSolver, typename TSharedRZSolver> template <typename TSolverIt> void SparseCoordinatorBase<TSharedZSolver,TSharedRZSolver>::setupStorage(Equations::SharedIEquation spEq, const SpectralFieldId id, const TSolverIt solveIt)
   {
      // Number of linear systems
      int nSystems = spEq->couplingInfo(id.second).nSystems();

      // Start row for storage information
      ArrayI startRow(nSystems);

      // Initialise the linear solver
      if((*solveIt)->nSystem() == 0)
      {
         // Initialise field storage and information
         for(int i = 0; i < nSystems; i++)
         {
            // Create data storage
            (*solveIt)->addStorage(spEq->couplingInfo(id.second).systemN(i), spEq->couplingInfo(id.second).rhsCols(i));
         }
      }

      // Build the solver matrices
      for(int i = 0; i < nSystems; i++)
      {
         // Store the start row
         startRow(i) = spEq->couplingInfo(id.second).fieldIndex()*spEq->couplingInfo(id.second).blockN(i);
      }

      // Store storage information
      (*solveIt)->addInformation(id,startRow);
   }

   template <typename TSharedZSolver, typename TSharedRZSolver> void SparseCoordinatorBase<TSharedZSolver,TSharedRZSolver>::transferOutput(const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq)
   {
      // Storage for identity
      SpectralFieldId myId;

      // Loop over all scalar equations
      std::vector<Equations::SharedIScalarEquation>::const_iterator scalEqIt;
      for(scalEqIt = scalEq.first; scalEqIt != scalEq.second; scalEqIt++)
      {
         // Get field identity
         myId = std::make_pair((*scalEqIt)->name(), FieldComponents::Spectral::SCALAR);

         // Get index of solver
         int myIdx = (*scalEqIt)->couplingInfo(myId.second).solverIndex();

         // Linear solve matrices are complex
         if((*scalEqIt)->couplingInfo(myId.second).isComplex())
         {
            // Create iterator to current complex solver
            SolverZ_iterator solZIt = this->mZSolvers.begin();
            std::advance(solZIt, myIdx);

            // Get solver output
            for(int i = 0; i < (*solZIt)->nSystem(); i++)
            {
               (*scalEqIt)->storeSolution(myId.second, (*solZIt)->solution(i), i, (*solZIt)->startRow(myId,i));
            }

         // Linear solve matrices are real
         } else
         {
            // Create iterator to current real solver
            SolverRZ_iterator solRZIt = this->mRZSolvers.begin();
            std::advance(solRZIt, myIdx);

            // Get solver output
            for(int i = 0; i < (*solRZIt)->nSystem(); i++)
            {
               (*scalEqIt)->storeSolution(myId.second, (*solRZIt)->solution(i), i, (*solRZIt)->startRow(myId,i));
            }
         }
      }

      // Loop over all vector equations
      std::vector<Equations::SharedIVectorEquation>::const_iterator vectEqIt;
      for(vectEqIt = vectEq.first; vectEqIt != vectEq.second; vectEqIt++)
      {
         // Loop over the vector equation components
         Equations::IVectorEquation::SpectralComponent_iterator compIt;
         Equations::IVectorEquation::SpectralComponent_range  compRange = (*vectEqIt)->spectralRange();
         for(compIt = compRange.first; compIt != compRange.second; ++compIt)
         {
            // Get field identity
            myId = std::make_pair((*vectEqIt)->name(), *compIt);

            // Get index of solver
            int myIdx = (*vectEqIt)->couplingInfo(myId.second).solverIndex();

            // Linear solve matrices are complex
            if((*vectEqIt)->couplingInfo(myId.second).isComplex())
            {
               // Create iterator to current complex solver
               SolverZ_iterator solZIt = this->mZSolvers.begin();
               std::advance(solZIt, myIdx);

               // Get solver output for first component
               for(int i = 0; i < (*solZIt)->nSystem(); i++)
               {
                  (*vectEqIt)->storeSolution(myId.second, (*solZIt)->solution(i), i, (*solZIt)->startRow(myId,i));
               }

               // Linear solve matrices are real
            } else
            {
               // Create iterator to current real solver
               SolverRZ_iterator solRZIt = this->mRZSolvers.begin();
               std::advance(solRZIt, myIdx);

               // Get solver output for first component
               for(int i = 0; i < (*solRZIt)->nSystem(); i++)
               {
                  (*vectEqIt)->storeSolution(myId.second, (*solRZIt)->solution(i), i, (*solRZIt)->startRow(myId,i));
               }
            }
         }
      }
   }

   template <typename TSharedZSolver, typename TSharedRZSolver> void SparseCoordinatorBase<TSharedZSolver,TSharedRZSolver>::initSolution(const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq)
   {
      // Storage for information and identity
      SpectralFieldId myId;

      // Loop over all scalar equations
      std::vector<Equations::SharedIScalarEquation>::const_iterator scalEqIt;
      for(scalEqIt = scalEq.first; scalEqIt != scalEq.second; scalEqIt++)
      {
         // Get field identity
         myId = std::make_pair((*scalEqIt)->name(), FieldComponents::Spectral::SCALAR);

         // Get index of solver
         int myIdx = (*scalEqIt)->couplingInfo(myId.second).solverIndex();

         // Linear solve matrices are complex
         if((*scalEqIt)->couplingInfo(myId.second).isComplex())
         {
            // Create iterator to current complex solver
            SolverZ_iterator solZIt = this->mZSolvers.begin();
            std::advance(solZIt, myIdx);

            // Get solver input
            for(int i = 0; i < (*solZIt)->nSystem(); i++)
            {
               Equations::copyUnknown(*(*scalEqIt), myId.second, (*solZIt)->rSolution(i), i, (*solZIt)->startRow(myId,i));
            }

         // Linear solve matrices are real
         } else
         {
            // Create iterator to current real solver
            SolverRZ_iterator solRZIt = this->mRZSolvers.begin();
            std::advance(solRZIt, myIdx);

            // Get solver input
            for(int i = 0; i < (*solRZIt)->nSystem(); i++)
            {
               Equations::copyUnknown(*(*scalEqIt), myId.second, (*solRZIt)->rSolution(i), i, (*solRZIt)->startRow(myId,i));
            }
         }
      }

      // Loop over all vector equations
      std::vector<Equations::SharedIVectorEquation>::const_iterator vectEqIt;
      for(vectEqIt = vectEq.first; vectEqIt != vectEq.second; vectEqIt++)
      {
         Equations::IVectorEquation::SpectralComponent_iterator compIt;
         Equations::IVectorEquation::SpectralComponent_range  compRange = (*vectEqIt)->spectralRange();
         for(compIt = compRange.first; compIt != compRange.second; ++compIt)
         {
            // Get field identity
            myId = std::make_pair((*vectEqIt)->name(), *compIt);

            // Get index of solver
            int myIdx = (*vectEqIt)->couplingInfo(myId.second).solverIndex();

            // Linear solve matrices are complex
            if((*vectEqIt)->couplingInfo(myId.second).isComplex())
            {
               // Create iterator to current complex solver
               SolverZ_iterator solZIt = this->mZSolvers.begin();
               std::advance(solZIt, myIdx);

               // Get solver input for toroidal component
               for(int i = 0; i < (*solZIt)->nSystem(); i++)
               {
                  Equations::copyUnknown(*(*vectEqIt), myId.second, (*solZIt)->rSolution(i), i, (*solZIt)->startRow(myId,i));
               }

               // Linear solve matrices are real
            } else
            {
               // Create iterator to current real solver
               SolverRZ_iterator solRZIt = this->mRZSolvers.begin();
               std::advance(solRZIt, myIdx);

               // Get solver input for toroidal component
               for(int i = 0; i < (*solRZIt)->nSystem(); i++)
               {
                  Equations::copyUnknown(*(*vectEqIt), myId.second, (*solRZIt)->rSolution(i), i, (*solRZIt)->startRow(myId,i));
               }
            }
         }
      }
   }

   template <typename TSharedZSolver, typename TSharedRZSolver> void SparseCoordinatorBase<TSharedZSolver,TSharedRZSolver>::getInput(const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq, const ScalarVariable_map& scalVar, const VectorVariable_map& vectVar)
   {
      // Storage for information and identity
      SpectralFieldId myId;

      // Loop over all scalar equations
      std::vector<Equations::SharedIScalarEquation>::const_iterator scalEqIt;
      for(scalEqIt = scalEq.first; scalEqIt != scalEq.second; scalEqIt++)
      {
         // Get field identity
         myId = std::make_pair((*scalEqIt)->name(), FieldComponents::Spectral::SCALAR);

         // Get index of solver
         int myIdx = (*scalEqIt)->couplingInfo(myId.second).solverIndex();

         // Linear solve matrices are complex
         if((*scalEqIt)->couplingInfo(myId.second).isComplex())
         {
            // Create iterator to current complex solver
            SolverZ_iterator solZIt = this->mZSolvers.begin();
            std::advance(solZIt, myIdx);

            // Get solver input
            this->getSolverInput(scalEqIt, myId, solZIt, scalVar, vectVar);

         // Linear solve matrices are real
         } else
         {
            // Create iterator to current real solver
            SolverRZ_iterator solRZIt = this->mRZSolvers.begin();
            std::advance(solRZIt, myIdx);

            // Get solver input
            this->getSolverInput(scalEqIt, myId, solRZIt, scalVar, vectVar);
         }
      }

      // Loop over all vector equations
      std::vector<Equations::SharedIVectorEquation>::const_iterator vectEqIt;
      for(vectEqIt = vectEq.first; vectEqIt != vectEq.second; vectEqIt++)
      {
         Equations::IVectorEquation::SpectralComponent_iterator compIt;
         Equations::IVectorEquation::SpectralComponent_range  compRange = (*vectEqIt)->spectralRange();
         for(compIt = compRange.first; compIt != compRange.second; ++compIt)
         {
            // Get field identity for first component
            myId = std::make_pair((*vectEqIt)->name(), *compIt);

            // Get index of solver
            int myIdx = (*vectEqIt)->couplingInfo(myId.second).solverIndex();

            // Linear solve matrices are complex
            if((*vectEqIt)->couplingInfo(myId.second).isComplex())
            {
               // Create iterator to current complex solver
               SolverZ_iterator solZIt = this->mZSolvers.begin();
               std::advance(solZIt, myIdx);

               // Get solver input
               this->getSolverInput(vectEqIt, myId, solZIt, scalVar, vectVar);

               // Linear solve matrices are real
            } else
            {
               // Create iterator to current real solver
               SolverRZ_iterator solRZIt = this->mRZSolvers.begin();
               std::advance(solRZIt, myIdx);

               // Get solver input
               this->getSolverInput(vectEqIt, myId, solRZIt, scalVar, vectVar);
            }
         }
      }
   }

   template <typename TSharedZSolver, typename TSharedRZSolver> template <typename TEquationIt, typename TSolverIt> void SparseCoordinatorBase<TSharedZSolver,TSharedRZSolver>::getSolverInput(const TEquationIt eqIt, const SpectralFieldId id, const TSolverIt solveIt, const ScalarVariable_map& scalVar, const VectorVariable_map& vectVar)
   {
      // Get timestep input
      for(int i = 0; i < (*solveIt)->nSystem(); i++)
      {
         // Copy field values into solver input
         Equations::copyNonlinear(*(*eqIt), id.second, (*solveIt)->rRHSData(i), i, (*solveIt)->startRow(id,i));

         // Apply quasi-inverse to nonlinear terms
         Equations::applyQuasiInverse(*(*eqIt), id.second, (*solveIt)->rRHSData(i), i, (*solveIt)->startRow(id,i));

         // Apply quasi-inverse to nonlinear terms
         Equations::addSource(*(*eqIt), id.second, (*solveIt)->rRHSData(i), i, (*solveIt)->startRow(id,i));

         // Loop over explicit fields
         Equations::CouplingInformation::FieldId_range   fRange = (*eqIt)->couplingInfo(id.second).explicitRange();
         for(Equations::CouplingInformation::FieldId_iterator  fIt = fRange.first; fIt != fRange.second; ++fIt)
         {
            if(fIt->second == FieldComponents::Spectral::SCALAR)
            {
               Equations::addExplicitLinear(*(*eqIt), id.second, (*solveIt)->rRHSData(i), (*solveIt)->startRow(id,i), *fIt, scalVar.find(fIt->first)->second->dom(0).perturbation(), i);
            } else
            {
               Equations::addExplicitLinear(*(*eqIt), id.second, (*solveIt)->rRHSData(i), (*solveIt)->startRow(id,i), *fIt, vectVar.find(fIt->first)->second->dom(0).perturbation().comp(fIt->second), i);
            }
         }
      }
   }
}
}

#endif // SPARSECOORDINATORBASE_HPP
