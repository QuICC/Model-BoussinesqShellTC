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

#include <iostream>

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
         
      protected:
         /**
          * @brief Create a real operator
          */
         void addSolver(RealSolver_iterator solIt, const int idx, const int start);

         /**
          * @brief Create a complex operator
          */
         void addSolver(ComplexSolver_iterator solIt, const int idx, const int start);

         /**
          * @brief Vector of (coupled) real operator
          */
         std::vector<SharedRealSolverType> mRealSolvers;

         /**
          * @brief Vector of (coupled) complex operator
          */
         std::vector<SharedComplexSolverType> mComplexSolvers;

      private:
   };

   /**
    * @brief Generic implementation to setup the solver storage
    */
   template <template <class,class> class TSol,typename TSolverIt> void setupSolverStorage(SparseCoordinatorData<TSol>& coord, Equations::SharedIEquation spEq, const int idx, SpectralFieldId);

   /**
    * @brief Generic implementation to store the solver solution
    */
   template <template <class,class> class TSol,typename TSolverIt, typename TEquation> void storeSolverSolution(SparseCoordinatorData<TSol>& coord, TEquation spEq, const int idx, SpectralFieldId);

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
    * @brief Generic implementation to compute RHS solvers
    */
   template <template <class,class> class TSol,typename TSolverIt> void computeRHSSolvers(SparseCoordinatorData<TSol>& coord, const int step);

   /**
    * @brief Generic implementation to update the time matrix solvers
    */
   template <template <class,class> class TSol,typename TSolverIt> void updateTimeMatrixSolvers(SparseCoordinatorData<TSol>& coord, const MHDFloat lhsCoeff, const MHDFloat rhsCoeff, const int step);

   /**
    * @brief Generic implementation to initialise the solver solution
    */
   template <template <class,class> class TSol,typename TSolverIt,typename TEquation> void initSolverSolution(SparseCoordinatorData<TSol>& coord, TEquation spEq, const int idx, SpectralFieldId id);

   /**
    * @brief Generic implementation to get the solver input
    */
   template <template <class,class> class TSol,typename TSolverIt,typename TEquation> void getSolverInput(SparseCoordinatorData<TSol>& coord, TEquation spEq, const int idx, SpectralFieldId id, const typename SparseCoordinatorData<TSol>::ScalarVariable_map& scalVar, const typename SparseCoordinatorData<TSol>::VectorVariable_map& vectVar);

   /**
    * @brief Compute the solver input independently of solver type
    */
   template <typename TEquationIt, typename TSolverIt> void computeSolverInput(const TEquationIt eqIt, const SpectralFieldId id, const TSolverIt solveIt, const std::map<PhysicalNames::Id, Datatypes::SharedScalarVariableType>& scalVar, const std::map<PhysicalNames::Id, Datatypes::SharedVectorVariableType>& vectVar);

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

   template <template <class,class> class TSolver> void SparseCoordinatorData<TSolver>::addSolver(typename SparseCoordinatorData<TSolver>::RealSolver_iterator solIt, const int idx, const int start)
   {
      if(idx > static_cast<int>(this->mRealSolvers.size()) - 1)
      {
         SparseCoordinatorData<TSolver>::SharedRealSolverType spSolver(new SparseCoordinatorData<TSolver>::RealSolverType(start));

         this->mRealSolvers.push_back(spSolver);
      }
   }

   template <template <class,class> class TSolver> void SparseCoordinatorData<TSolver>::addSolver(typename SparseCoordinatorData<TSolver>::ComplexSolver_iterator solIt, const int idx, const int start)
   {
      if(idx > static_cast<int>(this->mComplexSolvers.size()) - 1)
      {
         SparseCoordinatorData<TSolver>::SharedComplexSolverType spSolver(new SparseCoordinatorData<TSolver>::ComplexSolverType(start));

         this->mComplexSolvers.push_back(spSolver);
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

   template <template <class,class> class TSolver,typename TSolverIt,typename TEquation> void storeSolverSolution(SparseCoordinatorData<TSolver>& coord, TEquation spEq, const int idx, SpectralFieldId id)
   {
      // Create iterator to current solver
      TSolverIt solIt;
      coord.setIterator(solIt);
      std::advance(solIt, idx);

      // Get solver output
      for(int i = 0; i < (*solIt)->nSystem(); i++)
      {
         (*spEq)->storeSolution(id.second, (*solIt)->solution(i), i, (*solIt)->startRow(id,i));
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
         // Compute linear solve RHS
         (*solIt)->solve(step);
      }
   }

   template <template <class,class> class TSolver,typename TSolverIt> void computeRHSSolvers(SparseCoordinatorData<TSolver>& coord, const int step)
   {
      // Create iterator to current solver
      TSolverIt solIt;
      coord.setIterator(solIt);
      TSolverIt endIt;
      coord.setEndIterator(endIt);

      for(; solIt != endIt; ++solIt)
      {
         // Compute linear solve RHS
         (*solIt)->computeRHS(step);
      }
   }

   template <template <class,class> class TSolver,typename TSolverIt> void updateTimeMatrixSolvers(SparseCoordinatorData<TSolver>& coord, const MHDFloat lhsCoeff, const MHDFloat rhsCoeff, const int step)
   {
      // Create iterator to current solver
      TSolverIt solIt;
      coord.setIterator(solIt);
      TSolverIt endIt;
      coord.setEndIterator(endIt);

      for(; solIt != endIt; ++solIt)
      {
         // Compute linear solve RHS
         (*solIt)->updateTimeMatrix(lhsCoeff, rhsCoeff, step);
      }
   }

   template <template <class,class> class TSolver,typename TSolverIt,typename TEquation> void initSolverSolution(SparseCoordinatorData<TSolver>& coord, TEquation spEq, const int idx, SpectralFieldId id)
   {
      // Create iterator to current solver
      TSolverIt solIt;
      coord.setIterator(solIt);
      std::advance(solIt, idx);

      // Get solver input
      for(int i = 0; i < (*solIt)->nSystem(); i++)
      {
         Equations::copyUnknown(*(*spEq), id.second, (*solIt)->rSolution(i), i, (*solIt)->startRow(id,i));
      }
   }

   template <template <class,class> class TSolver,typename TSolverIt,typename TEquation> void getSolverInput(SparseCoordinatorData<TSolver>& coord, TEquation spEq, const int idx, SpectralFieldId id, const typename SparseCoordinatorData<TSolver>::ScalarVariable_map& scalVar, const typename SparseCoordinatorData<TSolver>::VectorVariable_map& vectVar)
   {
      // Create iterator to current complex field solver
      TSolverIt solIt;
      coord.setIterator(solIt);
      std::advance(solIt, idx);

      // Get solver input
      computeSolverInput(spEq, id, solIt, scalVar, vectVar);
   }

   //
   //
   //

   template <typename TEquationIt, typename TSolverIt> void computeSolverInput(const TEquationIt eqIt, const SpectralFieldId id, const TSolverIt solveIt, const std::map<PhysicalNames::Id, Datatypes::SharedScalarVariableType>& scalVar, const std::map<PhysicalNames::Id, Datatypes::SharedVectorVariableType>& vectVar)
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
         startRow(i) = spEq->couplingInfo(id.second).fieldIndex()*spEq->couplingInfo(id.second).blockN(i);
      }

      // Store storage information
      (*solveIt)->addInformation(id,startRow,spEq->pyName());
   }

}
}

#endif // SPARSECOORDINATORDATA_HPP
