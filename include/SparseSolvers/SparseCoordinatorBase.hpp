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
#include "Base/Typedefs.hpp"
#include "TypeSelectors/ScalarSelector.hpp"
#include "Equations/IScalarEquation.hpp"
#include "Equations/IVectorEquation.hpp"
#include "SparseSolvers/SparseCoordinatorData.hpp"

namespace GeoMHDiSCC {

namespace Solver {

   /**
    * @brief Implementation of the base for a general sparse solver coordinator
    */
   template <template <class,class> class TSolver> class SparseCoordinatorBase: public SparseCoordinatorData<TSolver>
   {
      public:
         /// Typedef for a shared scalar equation iterator
         typedef std::vector<Equations::SharedIScalarEquation>::iterator   ScalarEquation_iterator;

         /// Typedef for a shared vector equation iterator
         typedef std::vector<Equations::SharedIVectorEquation>::iterator   VectorEquation_iterator;

         /// Typedef for a shared scalar equation range
         typedef std::pair<ScalarEquation_iterator, ScalarEquation_iterator>  ScalarEquation_range;

         /// Typedef for a shared vector equation range
         typedef std::pair<VectorEquation_iterator, VectorEquation_iterator>  VectorEquation_range;

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
          * @brief Create the correct solver
          */
         void createSolver(Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp);

         /**
          * @brief Create storage for the solvers
          */
         void createStorage(Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp);

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
         void getInput(const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq, const typename SparseCoordinatorBase<TSolver>::ScalarVariable_map& scalVar, const typename SparseCoordinatorBase<TSolver>::VectorVariable_map& vectVar);

         /**
          * @brief Update equation unkowns with solver output 
          *
          * @param scalEq Shared scalar equations
          * @param vectEq Shared vector equations
          */
         void transferOutput(const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq);

         /**
          * @brief Update the step, counting from 0 to Nstep - 1
          */
         void updateStep();

         /**
          * @brief Number of solver substeps
          */
         int   mNStep;

         /**
          * @brief Current solver step
          */
         int   mStep;

      private:
   };

   template <template <class,class> class TSolver> SparseCoordinatorBase<TSolver>::SparseCoordinatorBase()
      : SparseCoordinatorData<TSolver>(), mNStep(1), mStep(0)
   {
   }

   template <template <class,class> class TSolver> SparseCoordinatorBase<TSolver>::~SparseCoordinatorBase()
   {
   }

   template <template <class,class> class TSolver> bool SparseCoordinatorBase<TSolver>::finishedStep() const
   {
      return (this->mStep == 0);
   }

   template <template <class,class> class TSolver> void SparseCoordinatorBase<TSolver>::updateStep()
   {
      if(this->solveTime() == SolveTiming::PROGNOSTIC || this->solveTime() == SolveTiming::AFTER)
      {
         this->mStep = (this->mStep + 1) % this->mNStep;
      }
   }

   template <template <class,class> class TSolver> void SparseCoordinatorBase<TSolver>::createSolver(Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp)
   {
      // System has a complex operator
      if(spEq->couplingInfo(comp).isComplex())
      {
         typename SparseCoordinatorBase<TSolver>::ComplexSolver_iterator solIt;
         this->addSolver(solIt, spEq->couplingInfo(comp).solverIndex(), spEq->couplingInfo(comp).fieldStart(), spEq->solveTiming());

      // System has a real operator
      } else
      {
         typename SparseCoordinatorBase<TSolver>::RealSolver_iterator solIt;
         this->addSolver(solIt, spEq->couplingInfo(comp).solverIndex(), spEq->couplingInfo(comp).fieldStart(), spEq->solveTiming());
      }
   }

   template <template <class,class> class TSolver> void SparseCoordinatorBase<TSolver>::createStorage(Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp)
   {
      // ID of the current field
      SpectralFieldId myId = std::make_pair(spEq->name(),comp);

      // Index of solver
      int myIdx = spEq->couplingInfo(myId.second).solverIndex();

      // System has a complex operator
      if(spEq->couplingInfo(myId.second).isComplex())
      {
         setupSolverStorage<TSolver,typename SparseCoordinatorBase<TSolver>::ComplexSolver_iterator>(*this, spEq, myIdx, myId);

      // System has a real operator
      } else
      {
         setupSolverStorage<TSolver,typename SparseCoordinatorBase<TSolver>::RealSolver_iterator>(*this, spEq, myIdx, myId);
      }
   }

   template <template <class,class> class TSolver> void SparseCoordinatorBase<TSolver>::transferOutput(const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq)
   {
      // Storage for identity
      SpectralFieldId myId;

      // Loop over all scalar equations
      std::vector<Equations::SharedIScalarEquation>::const_iterator scalEqIt;
      for(scalEqIt = scalEq.first; scalEqIt != scalEq.second; scalEqIt++)
      {
         if((*scalEqIt)->solveTiming() == this->solveTime())
         {
            // Get field identity
            myId = std::make_pair((*scalEqIt)->name(), FieldComponents::Spectral::SCALAR);

            // Get index of solver
            int myIdx = (*scalEqIt)->couplingInfo(myId.second).solverIndex();

            // System operator is complex
            if((*scalEqIt)->couplingInfo(myId.second).isComplex())
            {
               storeSolverSolution<TSolver,typename SparseCoordinatorBase<TSolver>::ComplexSolver_iterator>(*this, scalEqIt, myIdx, myId);

               // System operator is real real
            } else
            {
               storeSolverSolution<TSolver,typename SparseCoordinatorBase<TSolver>::RealSolver_iterator>(*this, scalEqIt, myIdx, myId);
            }
         }
      }

      // Loop over all vector equations
      std::vector<Equations::SharedIVectorEquation>::const_iterator vectEqIt;
      for(vectEqIt = vectEq.first; vectEqIt != vectEq.second; vectEqIt++)
      {
         if((*vectEqIt)->solveTiming() == this->solveTime())
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

               // System operator is complex
               if((*vectEqIt)->couplingInfo(myId.second).isComplex())
               {
                  storeSolverSolution<TSolver,typename SparseCoordinatorBase<TSolver>::ComplexSolver_iterator>(*this, vectEqIt, myIdx, myId);

                  // System operator is real
               } else
               {
                  storeSolverSolution<TSolver,typename SparseCoordinatorBase<TSolver>::RealSolver_iterator>(*this, vectEqIt, myIdx, myId);
               }
            }
         }
      }
   }

   template <template <class,class> class TSolver> void SparseCoordinatorBase<TSolver>::initSolution(const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq)
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

         // System operator is complex
         if((*scalEqIt)->couplingInfo(myId.second).isComplex())
         {
            initSolverSolution<TSolver,typename SparseCoordinatorBase<TSolver>::ComplexSolver_iterator>(*this, scalEqIt, myIdx, myId);

         // System operator is real
         } else
         {
            initSolverSolution<TSolver,typename SparseCoordinatorBase<TSolver>::RealSolver_iterator>(*this, scalEqIt, myIdx, myId);
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

            // System operator is complex
            if((*vectEqIt)->couplingInfo(myId.second).isComplex())
            {
               initSolverSolution<TSolver,typename SparseCoordinatorBase<TSolver>::ComplexSolver_iterator>(*this, vectEqIt, myIdx, myId);

            // System operator is real
            } else
            {
               initSolverSolution<TSolver,typename SparseCoordinatorBase<TSolver>::RealSolver_iterator>(*this, vectEqIt, myIdx, myId);
            }
         }
      }
   }

   template <template <class,class> class TSolver> void SparseCoordinatorBase<TSolver>::getInput(const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq, const typename SparseCoordinatorBase<TSolver>::ScalarVariable_map& scalVar, const typename SparseCoordinatorBase<TSolver>::VectorVariable_map& vectVar)
   {
      // Storage for information and identity
      SpectralFieldId myId;

      // Loop over all scalar equations
      std::vector<Equations::SharedIScalarEquation>::const_iterator scalEqIt;
      for(scalEqIt = scalEq.first; scalEqIt != scalEq.second; scalEqIt++)
      {
         if((*scalEqIt)->solveTiming() == this->solveTime())
         {
            // Get field identity
            myId = std::make_pair((*scalEqIt)->name(), FieldComponents::Spectral::SCALAR);

            // Get index of solver
            int myIdx = (*scalEqIt)->couplingInfo(myId.second).solverIndex();

            // System operator is complex
            if((*scalEqIt)->couplingInfo(myId.second).isComplex())
            {
               getSolverInput<TSolver,typename SparseCoordinatorBase<TSolver>::ComplexSolver_iterator>(*this, scalEqIt, myIdx, myId, scalVar, vectVar);

               // System operator is real
            } else
            {
               getSolverInput<TSolver,typename SparseCoordinatorBase<TSolver>::RealSolver_iterator>(*this, scalEqIt, myIdx, myId, scalVar, vectVar);
            }
         }
      }

      // Loop over all vector equations
      std::vector<Equations::SharedIVectorEquation>::const_iterator vectEqIt;
      for(vectEqIt = vectEq.first; vectEqIt != vectEq.second; vectEqIt++)
      {
         if((*vectEqIt)->solveTiming() == this->solveTime())
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
                  getSolverInput<TSolver,typename SparseCoordinatorBase<TSolver>::ComplexSolver_iterator>(*this, vectEqIt, myIdx, myId, scalVar, vectVar);

                  // Linear solve matrices are real
               } else
               {
                  getSolverInput<TSolver,typename SparseCoordinatorBase<TSolver>::RealSolver_iterator>(*this, vectEqIt, myIdx, myId, scalVar, vectVar);
               }
            }
         }
      }
   }
}
}

#endif // SPARSECOORDINATORBASE_HPP
