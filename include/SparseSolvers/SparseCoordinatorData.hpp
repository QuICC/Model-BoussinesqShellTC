/** 
 * @file SparseCoordinatorData.hpp
 * @brief Implementation of the base for a general sparse solver coordinator
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef SPARSECOORDINATORDATA_HPP
#define SPARSECOORDINATORDATA_HPP

// Configuration includes
//
#include "SmartPointers/SharedPtrMacro.h"

// System includes
//

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "TypeSelectors/SparseSolverDataSelector.hpp"
#include "Equations/IScalarEquation.hpp"
#include "Equations/IVectorEquation.hpp"
#include "SparseSolvers/SparseLinearSolver.hpp"
#include "SparseSolvers/SparseTrivialSolver.hpp"
#include "Timesteppers/SparseTimestepper.hpp"

namespace GeoMHDiSCC {

namespace Solver {

   /**
    * @brief Implementation of the base for a general sparse solver coordinator
    */
   template <template <class,class> class TSolver> class SparseCoordinatorData
   {
      public:
         /// Typedef for a real operator solver
         typedef typename SparseSolverDataSelector<TSolver>::RealSolverType RealSolverType;

         /// Typedef for a complex operator solver
         typedef typename SparseSolverDataSelector<TSolver>::ComplexSolverType ComplexSolverType;

         /// Typedef for a shared real operator solver
         typedef typename SharedPtrMacro<RealSolverType >  SharedRealSolverType;

         /// Typedef for a shared complex operator solver
         typedef typename SharedPtrMacro<ComplexSolverType >  SharedComplexSolverType;

         /// Typedef for an iterator to a real operator solver
         typedef typename std::vector<SharedRealSolverType>::iterator   RealSolver_iterator;

         /// Typedef for an iterator to a complex operator solver
         typedef typename std::vector<SharedComplexSolverType>::iterator   ComplexSolver_iterator;

         /// Typedef for a shared scalar variable map
         typedef std::map<PhysicalNames::Id, Datatypes::SharedScalarVariableType>  ScalarVariable_map;

         /// Typedef for a shared vector variable map
         typedef std::map<PhysicalNames::Id, Datatypes::SharedVectorVariableType>  VectorVariable_map;

         /**
          * @brief Constructor
          */
         SparseCoordinatorData();

         /**
          * @brief Destructor
          */
         virtual ~SparseCoordinatorData();

         /**
          * @brief Set iterator for real operator solver
          *
          * \mhdBug Should not be public
          */
         void setIterator(RealSolver_iterator& solIt);

         /**
          * @brief Set iterator for complex operator solver
          *
          * \mhdBug Should not be public
          */
         void setIterator(ComplexSolver_iterator& solIt);

         /**
          * @brief Set end iterator for real operator solver
          *
          * \mhdBug Should not be public
          */
         void setEndIterator(RealSolver_iterator& solIt);

         /**
          * @brief Set end iterator for complex operator solver
          *
          * \mhdBug Should not be public
          */
         void setEndIterator(ComplexSolver_iterator& solIt);

         /**
          * @brief Get current solver time
          */
         SolveTiming::Id solveTime() const;

         /**
          * @brief Set solve time
          */
         void setSolveTime(const SolveTiming::Id time);

         /**
          * @brief Clear the RHS data of all solvers
          */
         void clearSolvers();
         
      protected:
         /**
          * @brief Create a real operator
          */
         void addSolver(RealSolver_iterator solIt, const int idx, const int start, const SolveTiming::Id time);

         /**
          * @brief Create a complex operator
          */
         void addSolver(ComplexSolver_iterator solIt, const int idx, const int start, const SolveTiming::Id time);

         /**
          * @brief Initi the start rows for the solvers
          */
         void initStartRow();

         /**
          * @brief Vector of (coupled) real operator
          */
         std::vector<SharedRealSolverType> mRealSolvers;

         /**
          * @brief Vector of (coupled) complex operator
          */
         std::vector<SharedComplexSolverType> mComplexSolvers;

         /**
          * @brief Storage for the current solve time
          */
         SolveTiming::Id   mSolveTime;

      private:
   };

   /**
    * @brief Generic implementation to setup the solver storage
    */
   template <template <class,class> class TSol,typename TSolverIt> void setupSolverStorage(SparseCoordinatorData<TSol>& coord, Equations::SharedIEquation spEq, const int idx, SpectralFieldId);

   /**
    * @brief Generic implementation to store the solver solution
    */
   template <template <class,class> class TSol,typename TSolverIt, typename TEquationIt> void storeSolverSolution(SparseCoordinatorData<TSol>& coord, TEquationIt eqIt, const int idx, SpectralFieldId);

   /**
    * @brief Generic implementation to initialise the solvers
    */
   template <template <class,class> class TSol,typename TSolverIt> void initSolvers(SparseCoordinatorData<TSol>& coord);

   /**
    * @brief Generic implementation to update the solvers
    */
   template <template <class,class> class TSol,typename TSolverIt> void updateSolvers(SparseCoordinatorData<TSol>& coord);

   /**
    * @brief Generic implementation to solver solvers
    */
   template <template <class,class> class TSol,typename TSolverIt> void solveSolvers(SparseCoordinatorData<TSol>& coord, const int step);

   /**
    * @brief Generic implementation to update the time matrix solvers
    */
   template <template <class,class> class TSol,typename TSolverIt> void updateTimeMatrixSolvers(SparseCoordinatorData<TSol>& coord, const int step, const MHDFloat dt);

   /**
    * @brief Generic implementation to initialise the solver solution
    */
   template <template <class,class> class TSol,typename TSolverIt,typename TEquationIt> void initSolverSolution(SparseCoordinatorData<TSol>& coord, TEquationIt eqIt, const int idx, SpectralFieldId id);

   /**
    * @brief Generic implementation to get the explicit linear input
    */
   template <template <class,class> class TSol,typename TSolverIt,typename TEquationIt> void getExplicitSolverInput(SparseCoordinatorData<TSol>& coord, TEquationIt eqIt, const int idx, SpectralFieldId id, const ModelOperator::Id opId, const typename SparseCoordinatorData<TSol>::ScalarVariable_map& scalVar, const typename SparseCoordinatorData<TSol>::VectorVariable_map& vectVar);

   /**
    * @brief Generic implementation to get the solver input
    */
   template <template <class,class> class TSol,typename TSolverIt,typename TEquationIt> void getSolverInput(SparseCoordinatorData<TSol>& coord, TEquationIt eqIt, const int idx, SpectralFieldId id, const typename SparseCoordinatorData<TSol>::ScalarVariable_map& scalVar, const typename SparseCoordinatorData<TSol>::VectorVariable_map& vectVar);

   /**
    * @brief Compute the explicit linear input independently of solver type
    */
   template <typename TSolverIt,typename TEquationIt> void computeExplicitSolverInput(const TEquationIt eqIt, const SpectralFieldId id, const TSolverIt solveIt, const ModelOperator::Id opId, const std::map<PhysicalNames::Id, Datatypes::SharedScalarVariableType>& scalVar, const std::map<PhysicalNames::Id, Datatypes::SharedVectorVariableType>& vectVar);

   /**
    * @brief Compute the solver input independently of solver type
    */
   template <typename TSolverIt,typename TEquationIt> void computeSolverInput(const TEquationIt eqIt, const SpectralFieldId id, const TSolverIt solveIt, const std::map<PhysicalNames::Id, Datatypes::SharedScalarVariableType>& scalVar, const std::map<PhysicalNames::Id, Datatypes::SharedVectorVariableType>& vectVar);

   /**
    * @brief Setup the storage and related information independently of solver type
    */
   template <typename TSolverIt> void setupStorage(Equations::SharedIEquation spEq, const SpectralFieldId id, const TSolverIt solveIt);

   //
   //
   //

   template <template <class,class> class TSolver> SparseCoordinatorData<TSolver>::SparseCoordinatorData()
   {
   }

   template <template <class,class> class TSolver> SparseCoordinatorData<TSolver>::~SparseCoordinatorData()
   {
   }

   template <template <class,class> class TSolver> SolveTiming::Id SparseCoordinatorData<TSolver>::solveTime() const
   {
      return this->mSolveTime;
   }

   template <template <class,class> class TSolver> void SparseCoordinatorData<TSolver>::setSolveTime(const SolveTiming::Id time)
   {
      this->mSolveTime = time;
   }

   template <template <class,class> class TSolver> void SparseCoordinatorData<TSolver>::addSolver(typename SparseCoordinatorData<TSolver>::RealSolver_iterator solIt, const int idx, const int start, const SolveTiming::Id time)
   {
      if(idx > static_cast<int>(this->mRealSolvers.size()) - 1)
      {
         SparseCoordinatorData<TSolver>::SharedRealSolverType spSolver(new SparseCoordinatorData<TSolver>::RealSolverType(start,time));

         this->mRealSolvers.push_back(spSolver);
      }
   }

   template <template <class,class> class TSolver> void SparseCoordinatorData<TSolver>::addSolver(typename SparseCoordinatorData<TSolver>::ComplexSolver_iterator solIt, const int idx, const int start, const SolveTiming::Id time)
   {
      if(idx > static_cast<int>(this->mComplexSolvers.size()) - 1)
      {
         SparseCoordinatorData<TSolver>::SharedComplexSolverType spSolver(new SparseCoordinatorData<TSolver>::ComplexSolverType(start,time));

         this->mComplexSolvers.push_back(spSolver);
      }
   }

   template <template <class,class> class TSolver> void SparseCoordinatorData<TSolver>::initStartRow()
   {
      for(typename std::vector<SharedRealSolverType>::iterator rIt = this->mRealSolvers.begin(); rIt != this->mRealSolvers.end(); ++rIt)
      {
         (*rIt)->initStartRow();
      }

      for(typename std::vector<SharedComplexSolverType>::iterator zIt = this->mComplexSolvers.begin(); zIt != this->mComplexSolvers.end(); ++zIt)
      {
         (*zIt)->initStartRow();
      }
   }

   template <template <class,class> class TSolver> void SparseCoordinatorData<TSolver>::clearSolvers()
   {
      for(typename std::vector<SharedRealSolverType>::iterator rIt = this->mRealSolvers.begin(); rIt != this->mRealSolvers.end(); ++rIt)
      {
         (*rIt)->zeroSolver();
      }

      for(typename std::vector<SharedComplexSolverType>::iterator zIt = this->mComplexSolvers.begin(); zIt != this->mComplexSolvers.end(); ++zIt)
      {
         (*zIt)->zeroSolver();
      }
   }

   template <template <class,class> class TSolver> void SparseCoordinatorData<TSolver>::setIterator(typename SparseCoordinatorData<TSolver>::RealSolver_iterator& solIt)
   {
      solIt = this->mRealSolvers.begin();
   }

   template <template <class,class> class TSolver> void SparseCoordinatorData<TSolver>::setIterator(typename SparseCoordinatorData<TSolver>::ComplexSolver_iterator& solIt)
   {
      solIt = this->mComplexSolvers.begin();
   }

   template <template <class,class> class TSolver> void SparseCoordinatorData<TSolver>::setEndIterator(typename SparseCoordinatorData<TSolver>::RealSolver_iterator& solIt)
   {
      solIt = this->mRealSolvers.end();
   }

   template <template <class,class> class TSolver> void SparseCoordinatorData<TSolver>::setEndIterator(typename SparseCoordinatorData<TSolver>::ComplexSolver_iterator& solIt)
   {
      solIt = this->mComplexSolvers.end();
   }

   //
   //
   //
   //

   template <template <class,class> class TSolver, typename TSolverIt> void setupSolverStorage(SparseCoordinatorData<TSolver>& coord, Equations::SharedIEquation spEq, const int idx, SpectralFieldId id)
   {
      // Create iterator to current solver
      TSolverIt solIt;
      coord.setIterator(solIt);
      std::advance(solIt, idx);

      // setup storage and information
      setupStorage(spEq, id, solIt);
   }

   template <template <class,class> class TSolver,typename TSolverIt,typename TEquationIt> void storeSolverSolution(SparseCoordinatorData<TSolver>& coord, TEquationIt eqIt, const int idx, SpectralFieldId id)
   {
      // Create iterator to current solver
      TSolverIt solIt;
      coord.setIterator(solIt);
      std::advance(solIt, idx);

      // Get solver output
      for(int i = 0; i < (*solIt)->nSystem(); i++)
      {
         (*eqIt)->storeSolution(id.second, (*solIt)->solution(i), i, (*solIt)->startRow(id,i));
      }
   }

   template <template <class,class> class TSolver,typename TSolverIt> void initSolvers(SparseCoordinatorData<TSolver>& coord)
   {
      // Create iterator to current solver
      TSolverIt solIt;
      coord.setIterator(solIt);
      TSolverIt endIt;
      coord.setEndIterator(endIt);

      for(; solIt != endIt; ++solIt)
      {
         (*solIt)->initSolver();
      }
   }

   template <template <class,class> class TSolver,typename TSolverIt> void updateSolvers(SparseCoordinatorData<TSolver>& coord)
   {
      // Create iterator to current solver
      TSolverIt solIt;
      coord.setIterator(solIt);
      TSolverIt endIt;
      coord.setEndIterator(endIt);

      for(; solIt != endIt; ++solIt)
      {
         (*solIt)->updateSolver();
      }
   }

   template <template <class,class> class TSolver,typename TSolverIt> void solveSolvers(SparseCoordinatorData<TSolver>& coord, const int step)
   {
      // Create iterator to current solver
      TSolverIt solIt;
      coord.setIterator(solIt);
      TSolverIt endIt;
      coord.setEndIterator(endIt);

      for(; solIt != endIt; ++solIt)
      {
         if((*solIt)->solveTiming() == coord.solveTime())
         {
            // Prepare solve of linear system
            bool needSolve = (*solIt)->solve(step);

            if(needSolve)
            {
               // Solve linear system
               (*solIt)->solve(step);

               // Work on fields after solve
               (*solIt)->postSolve(step);
            }
         }
      }
   }

   template <template <class,class> class TSolver,typename TSolverIt> void updateTimeMatrixSolvers(SparseCoordinatorData<TSolver>& coord, const int step, const MHDFloat dt)
   {
      // Create iterator to current solver
      TSolverIt solIt;
      coord.setIterator(solIt);
      TSolverIt endIt;
      coord.setEndIterator(endIt);

      for(; solIt != endIt; ++solIt)
      {
         // Compute linear solve RHS
         (*solIt)->updateTimeMatrix(step, dt);
      }
   }

   template <template <class,class> class TSolver,typename TSolverIt,typename TEquationIt> void initSolverSolution(SparseCoordinatorData<TSolver>& coord, TEquationIt eqIt, const int idx, SpectralFieldId id)
   {
      // Create iterator to current solver
      TSolverIt solIt;
      coord.setIterator(solIt);
      std::advance(solIt, idx);

      // Get solver input
      for(int i = 0; i < (*solIt)->nSystem(); i++)
      {
         if((*eqIt)->couplingInfo(id.second).isGalerkin())
         {
            Equations::solveStencilUnknown(*(*eqIt), id.second, (*solIt)->rSolution(i), i, (*solIt)->startRow(id,i));

         } else
         {
            Equations::copyUnknown(*(*eqIt), id.second, (*solIt)->rSolution(i), i, (*solIt)->startRow(id,i), true, true);
         }
      }

      (*solIt)->initSolutions();
   }

   template <template <class,class> class TSolver,typename TSolverIt,typename TEquationIt> void getExplicitSolverInput(SparseCoordinatorData<TSolver>& coord, TEquationIt eqIt, const int idx, SpectralFieldId id, const ModelOperator::Id opId, const typename SparseCoordinatorData<TSolver>::ScalarVariable_map& scalVar, const typename SparseCoordinatorData<TSolver>::VectorVariable_map& vectVar)
   {
      // Create iterator to current complex field solver
      TSolverIt solIt;
      coord.setIterator(solIt);
      std::advance(solIt, idx);

      // Get solver input
      computeExplicitSolverInput(eqIt, id, solIt, opId, scalVar, vectVar);
   }

   template <template <class,class> class TSolver,typename TSolverIt,typename TEquationIt> void getSolverInput(SparseCoordinatorData<TSolver>& coord, TEquationIt eqIt, const int idx, SpectralFieldId id, const typename SparseCoordinatorData<TSolver>::ScalarVariable_map& scalVar, const typename SparseCoordinatorData<TSolver>::VectorVariable_map& vectVar)
   {
      // Create iterator to current complex field solver
      TSolverIt solIt;
      coord.setIterator(solIt);
      std::advance(solIt, idx);

      // Get solver input
      computeSolverInput(eqIt, id, solIt, scalVar, vectVar);
   }

   //
   //
   //

   template <typename TSolverIt,typename TEquationIt> void computeExplicitSolverInput(const TEquationIt eqIt, const SpectralFieldId id, const TSolverIt solveIt, const ModelOperator::Id opId, const std::map<PhysicalNames::Id, Datatypes::SharedScalarVariableType>& scalVar, const std::map<PhysicalNames::Id, Datatypes::SharedVectorVariableType>& vectVar)
   {
      // Get timestep input
      for(int i = 0; i < (*solveIt)->nSystem(); i++)
      {
         // Loop over explicit fields
         Equations::CouplingInformation::FieldId_range   fRange = (*eqIt)->couplingInfo(id.second).explicitRange(opId);
         for(Equations::CouplingInformation::FieldId_iterator  fIt = fRange.first; fIt != fRange.second; ++fIt)
         {
            if(fIt->second == FieldComponents::Spectral::SCALAR)
            {
               Equations::addExplicitTerm(*(*eqIt), opId, id.second, (*solveIt)->rRHSData(i), (*solveIt)->startRow(id,i), *fIt, scalVar.find(fIt->first)->second->dom(0).perturbation(), i);
            } else
            {
               Equations::addExplicitTerm(*(*eqIt), opId, id.second, (*solveIt)->rRHSData(i), (*solveIt)->startRow(id,i), *fIt, vectVar.find(fIt->first)->second->dom(0).perturbation().comp(fIt->second), i);
            }
         }
      }
   }

   template <typename TSolverIt,typename TEquationIt> void computeSolverInput(const TEquationIt eqIt, const SpectralFieldId id, const TSolverIt solveIt, const std::map<PhysicalNames::Id, Datatypes::SharedScalarVariableType>& scalVar, const std::map<PhysicalNames::Id, Datatypes::SharedVectorVariableType>& vectVar)
   {
      // Get timestep input
      for(int i = 0; i < (*solveIt)->nSystem(); i++)
      {
         // Copy field values into solver input
         Equations::copyNonlinear(*(*eqIt), id.second, (*solveIt)->rRHSData(i), i, (*solveIt)->startRow(id,i));

         // Apply quasi-inverse to nonlinear terms
         Equations::addSource(*(*eqIt), id.second, (*solveIt)->rRHSData(i), i, (*solveIt)->startRow(id,i));
      }
   }

   template <typename TSolverIt> void setupStorage(Equations::SharedIEquation spEq, const SpectralFieldId id, const TSolverIt solveIt)
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

      // Gather information about the solver matrices
      for(int i = 0; i < nSystems; i++)
      {
         // Store the start row
         startRow(i) = spEq->couplingInfo(id.second).galerkinN(i);
      }

      // Store storage information
      (*solveIt)->addInformation(id, spEq->couplingInfo(id.second).fieldIndex(), startRow);
   }

   //
   // Dummy solver specializations
   //


   // Solver independent

   template <> inline void setupStorage<ComplexDummy_iterator>(Equations::SharedIEquation spEq, const SpectralFieldId id, const ComplexDummy_iterator solveIt) {};

   template <> inline void computeSolverInput<ComplexDummy_iterator,std::vector<Equations::SharedIScalarEquation>::const_iterator>(const std::vector<Equations::SharedIScalarEquation>::const_iterator eqIt, const SpectralFieldId id, const ComplexDummy_iterator solveIt, const std::map<PhysicalNames::Id, Datatypes::SharedScalarVariableType>& scalVar, const std::map<PhysicalNames::Id, Datatypes::SharedVectorVariableType>& vectVar) {};
   template <> inline void computeSolverInput<ComplexDummy_iterator,std::vector<Equations::SharedIVectorEquation>::const_iterator>(const std::vector<Equations::SharedIVectorEquation>::const_iterator eqIt, const SpectralFieldId id, const ComplexDummy_iterator solveIt, const std::map<PhysicalNames::Id, Datatypes::SharedScalarVariableType>& scalVar, const std::map<PhysicalNames::Id, Datatypes::SharedVectorVariableType>& vectVar) {};

   template <> inline void computeExplicitSolverInput<ComplexDummy_iterator,std::vector<Equations::SharedIScalarEquation>::const_iterator>(const std::vector<Equations::SharedIScalarEquation>::const_iterator eqIt, const SpectralFieldId id, const ComplexDummy_iterator solveIt, const ModelOperator::Id opId, const std::map<PhysicalNames::Id, Datatypes::SharedScalarVariableType>& scalVar, const std::map<PhysicalNames::Id, Datatypes::SharedVectorVariableType>& vectVar) {};
   template <> inline void computeExplicitSolverInput<ComplexDummy_iterator,std::vector<Equations::SharedIVectorEquation>::const_iterator>(const std::vector<Equations::SharedIVectorEquation>::const_iterator eqIt, const SpectralFieldId id, const ComplexDummy_iterator solveIt, const ModelOperator::Id opId, const std::map<PhysicalNames::Id, Datatypes::SharedScalarVariableType>& scalVar, const std::map<PhysicalNames::Id, Datatypes::SharedVectorVariableType>& vectVar) {};

   // SparseLinearSolver

   template <> inline void setupSolverStorage<SparseLinearSolver,ComplexDummy_iterator>(SparseCoordinatorData<SparseLinearSolver>& coord, Equations::SharedIEquation spEq, const int idx, SpectralFieldId) {};

   template <> inline void initSolvers<SparseLinearSolver,ComplexDummy_iterator>(SparseCoordinatorData<SparseLinearSolver>& coord) {};

   template <> inline void updateSolvers<SparseLinearSolver,ComplexDummy_iterator>(SparseCoordinatorData<SparseLinearSolver>& coord) {};

   template <> inline void solveSolvers<SparseLinearSolver,ComplexDummy_iterator>(SparseCoordinatorData<SparseLinearSolver>& coord, const int step) {};

   template <> inline void updateTimeMatrixSolvers<SparseLinearSolver,ComplexDummy_iterator>(SparseCoordinatorData<SparseLinearSolver>& coord, const int step, const MHDFloat dt) {};

   template <> inline void storeSolverSolution<SparseLinearSolver,ComplexDummy_iterator>(SparseCoordinatorData<SparseLinearSolver>& coord, std::vector<Equations::SharedIScalarEquation>::const_iterator eqIt, const int idx, SpectralFieldId) {};
   template <> inline void storeSolverSolution<SparseLinearSolver,ComplexDummy_iterator>(SparseCoordinatorData<SparseLinearSolver>& coord, std::vector<Equations::SharedIVectorEquation>::const_iterator eqIt, const int idx, SpectralFieldId) {};

   template <> inline void initSolverSolution<SparseLinearSolver,ComplexDummy_iterator>(SparseCoordinatorData<SparseLinearSolver>& coord, std::vector<Equations::SharedIScalarEquation>::const_iterator eqIt, const int idx, SpectralFieldId id) {};
   template <> inline void initSolverSolution<SparseLinearSolver,ComplexDummy_iterator>(SparseCoordinatorData<SparseLinearSolver>& coord, std::vector<Equations::SharedIVectorEquation>::const_iterator eqIt, const int idx, SpectralFieldId id) {};

   template <> inline void getExplicitSolverInput<SparseLinearSolver,ComplexDummy_iterator>(SparseCoordinatorData<SparseLinearSolver>& coord, std::vector<Equations::SharedIScalarEquation>::const_iterator eqIt, const int idx, SpectralFieldId id, const ModelOperator::Id opId, const SparseCoordinatorData<SparseLinearSolver>::ScalarVariable_map& scalVar, const SparseCoordinatorData<SparseLinearSolver>::VectorVariable_map& vectVar) {};
   template <> inline void getExplicitSolverInput<SparseLinearSolver,ComplexDummy_iterator>(SparseCoordinatorData<SparseLinearSolver>& coord, std::vector<Equations::SharedIVectorEquation>::const_iterator eqIt, const int idx, SpectralFieldId id, const ModelOperator::Id opId, const SparseCoordinatorData<SparseLinearSolver>::ScalarVariable_map& scalVar, const SparseCoordinatorData<SparseLinearSolver>::VectorVariable_map& vectVar) {};

   template <> inline void getSolverInput<SparseLinearSolver,ComplexDummy_iterator>(SparseCoordinatorData<SparseLinearSolver>& coord, std::vector<Equations::SharedIScalarEquation>::const_iterator eqIt, const int idx, SpectralFieldId id, const SparseCoordinatorData<SparseLinearSolver>::ScalarVariable_map& scalVar, const SparseCoordinatorData<SparseLinearSolver>::VectorVariable_map& vectVar) {};
   template <> inline void getSolverInput<SparseLinearSolver,ComplexDummy_iterator>(SparseCoordinatorData<SparseLinearSolver>& coord, std::vector<Equations::SharedIVectorEquation>::const_iterator eqIt, const int idx, SpectralFieldId id, const SparseCoordinatorData<SparseLinearSolver>::ScalarVariable_map& scalVar, const SparseCoordinatorData<SparseLinearSolver>::VectorVariable_map& vectVar) {};

   // SparseTrivialSolver

   template <> inline void setupSolverStorage<SparseTrivialSolver,ComplexDummy_iterator>(SparseCoordinatorData<SparseTrivialSolver>& coord, Equations::SharedIEquation spEq, const int idx, SpectralFieldId) {};

   template <> inline void initSolvers<SparseTrivialSolver,ComplexDummy_iterator>(SparseCoordinatorData<SparseTrivialSolver>& coord) {};

   template <> inline void updateSolvers<SparseTrivialSolver,ComplexDummy_iterator>(SparseCoordinatorData<SparseTrivialSolver>& coord) {};

   template <> inline void solveSolvers<SparseTrivialSolver,ComplexDummy_iterator>(SparseCoordinatorData<SparseTrivialSolver>& coord, const int step) {};

   template <> inline void updateTimeMatrixSolvers<SparseTrivialSolver,ComplexDummy_iterator>(SparseCoordinatorData<SparseTrivialSolver>& coord, const int step, const MHDFloat dt) {};

   template <> inline void storeSolverSolution<SparseTrivialSolver,ComplexDummy_iterator>(SparseCoordinatorData<SparseTrivialSolver>& coord, std::vector<Equations::SharedIScalarEquation>::const_iterator eqIt, const int idx, SpectralFieldId) {};
   template <> inline void storeSolverSolution<SparseTrivialSolver,ComplexDummy_iterator>(SparseCoordinatorData<SparseTrivialSolver>& coord, std::vector<Equations::SharedIVectorEquation>::const_iterator eqIt, const int idx, SpectralFieldId) {};

   template <> inline void initSolverSolution<SparseTrivialSolver,ComplexDummy_iterator>(SparseCoordinatorData<SparseTrivialSolver>& coord, std::vector<Equations::SharedIScalarEquation>::const_iterator eqIt, const int idx, SpectralFieldId id) {};
   template <> inline void initSolverSolution<SparseTrivialSolver,ComplexDummy_iterator>(SparseCoordinatorData<SparseTrivialSolver>& coord, std::vector<Equations::SharedIVectorEquation>::const_iterator eqIt, const int idx, SpectralFieldId id) {};

   template <> inline void getExplicitSolverInput<SparseTrivialSolver,ComplexDummy_iterator>(SparseCoordinatorData<SparseTrivialSolver>& coord, std::vector<Equations::SharedIScalarEquation>::const_iterator eqIt, const int idx, SpectralFieldId id, const ModelOperator::Id opId, const SparseCoordinatorData<SparseTrivialSolver>::ScalarVariable_map& scalVar, const SparseCoordinatorData<SparseTrivialSolver>::VectorVariable_map& vectVar) {};
   template <> inline void getExplicitSolverInput<SparseTrivialSolver,ComplexDummy_iterator>(SparseCoordinatorData<SparseTrivialSolver>& coord, std::vector<Equations::SharedIVectorEquation>::const_iterator eqIt, const int idx, SpectralFieldId id, const ModelOperator::Id opId, const SparseCoordinatorData<SparseTrivialSolver>::ScalarVariable_map& scalVar, const SparseCoordinatorData<SparseTrivialSolver>::VectorVariable_map& vectVar) {};

   template <> inline void getSolverInput<SparseTrivialSolver,ComplexDummy_iterator>(SparseCoordinatorData<SparseTrivialSolver>& coord, std::vector<Equations::SharedIScalarEquation>::const_iterator eqIt, const int idx, SpectralFieldId id, const SparseCoordinatorData<SparseTrivialSolver>::ScalarVariable_map& scalVar, const SparseCoordinatorData<SparseTrivialSolver>::VectorVariable_map& vectVar) {};
   template <> inline void getSolverInput<SparseTrivialSolver,ComplexDummy_iterator>(SparseCoordinatorData<SparseTrivialSolver>& coord, std::vector<Equations::SharedIVectorEquation>::const_iterator eqIt, const int idx, SpectralFieldId id, const SparseCoordinatorData<SparseTrivialSolver>::ScalarVariable_map& scalVar, const SparseCoordinatorData<SparseTrivialSolver>::VectorVariable_map& vectVar) {};

   // SparseTimestepper

   template <> inline void setupSolverStorage<Timestep::SparseTimestepper,ComplexDummy_iterator>(SparseCoordinatorData<Timestep::SparseTimestepper>& coord, Equations::SharedIEquation spEq, const int idx, SpectralFieldId) {};

   template <> inline void initSolvers<Timestep::SparseTimestepper,ComplexDummy_iterator>(SparseCoordinatorData<Timestep::SparseTimestepper>& coord) {};

   template <> inline void updateSolvers<Timestep::SparseTimestepper,ComplexDummy_iterator>(SparseCoordinatorData<Timestep::SparseTimestepper>& coord) {};

   template <> inline void solveSolvers<Timestep::SparseTimestepper,ComplexDummy_iterator>(SparseCoordinatorData<Timestep::SparseTimestepper>& coord, const int step) {};

   template <> inline void updateTimeMatrixSolvers<Timestep::SparseTimestepper,ComplexDummy_iterator>(SparseCoordinatorData<Timestep::SparseTimestepper>& coord, const int step, const MHDFloat dt) {};

   template <> inline void storeSolverSolution<Timestep::SparseTimestepper,ComplexDummy_iterator>(SparseCoordinatorData<Timestep::SparseTimestepper>& coord, std::vector<Equations::SharedIScalarEquation>::const_iterator eqIt, const int idx, SpectralFieldId) {};
   template <> inline void storeSolverSolution<Timestep::SparseTimestepper,ComplexDummy_iterator>(SparseCoordinatorData<Timestep::SparseTimestepper>& coord, std::vector<Equations::SharedIVectorEquation>::const_iterator eqIt, const int idx, SpectralFieldId) {};

   template <> inline void initSolverSolution<Timestep::SparseTimestepper,ComplexDummy_iterator>(SparseCoordinatorData<Timestep::SparseTimestepper>& coord, std::vector<Equations::SharedIScalarEquation>::const_iterator eqIt, const int idx, SpectralFieldId id) {};
   template <> inline void initSolverSolution<Timestep::SparseTimestepper,ComplexDummy_iterator>(SparseCoordinatorData<Timestep::SparseTimestepper>& coord, std::vector<Equations::SharedIVectorEquation>::const_iterator eqIt, const int idx, SpectralFieldId id) {};

   template <> inline void getExplicitSolverInput<Timestep::SparseTimestepper,ComplexDummy_iterator>(SparseCoordinatorData<Timestep::SparseTimestepper>& coord, std::vector<Equations::SharedIScalarEquation>::const_iterator eqIt, const int idx, SpectralFieldId id, const ModelOperator::Id opId, const SparseCoordinatorData<Timestep::SparseTimestepper>::ScalarVariable_map& scalVar, const SparseCoordinatorData<Timestep::SparseTimestepper>::VectorVariable_map& vectVar) {};
   template <> inline void getExplicitSolverInput<Timestep::SparseTimestepper,ComplexDummy_iterator>(SparseCoordinatorData<Timestep::SparseTimestepper>& coord, std::vector<Equations::SharedIVectorEquation>::const_iterator eqIt, const int idx, SpectralFieldId id, const ModelOperator::Id opId, const SparseCoordinatorData<Timestep::SparseTimestepper>::ScalarVariable_map& scalVar, const SparseCoordinatorData<Timestep::SparseTimestepper>::VectorVariable_map& vectVar) {};

   template <> inline void getSolverInput<Timestep::SparseTimestepper,ComplexDummy_iterator>(SparseCoordinatorData<Timestep::SparseTimestepper>& coord, std::vector<Equations::SharedIScalarEquation>::const_iterator eqIt, const int idx, SpectralFieldId id, const SparseCoordinatorData<Timestep::SparseTimestepper>::ScalarVariable_map& scalVar, const SparseCoordinatorData<Timestep::SparseTimestepper>::VectorVariable_map& vectVar) {};
   template <> inline void getSolverInput<Timestep::SparseTimestepper,ComplexDummy_iterator>(SparseCoordinatorData<Timestep::SparseTimestepper>& coord, std::vector<Equations::SharedIVectorEquation>::const_iterator eqIt, const int idx, SpectralFieldId id, const SparseCoordinatorData<Timestep::SparseTimestepper>::ScalarVariable_map& scalVar, const SparseCoordinatorData<Timestep::SparseTimestepper>::VectorVariable_map& vectVar) {};


}
}

#endif // SPARSECOORDINATORDATA_HPP
