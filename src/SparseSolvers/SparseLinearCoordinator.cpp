/** \file SparseLinearCoordinator.cpp
 *  \brief Implementation of a general linear solver structure
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
#include "SparseSolvers/SparseLinearCoordinator.hpp"

// Project includes
//
#include "Base/MathConstants.hpp"

namespace GeoMHDiSCC {

namespace Solver {

   SparseLinearCoordinator::SparseLinearCoordinator()
      : SparseLinearCoordinatorBase()
   {
   }

   SparseLinearCoordinator::~SparseLinearCoordinator()
   {
   }

   void SparseLinearCoordinator::solve(const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq)
   {
      // Update the equation input to the timestepper
      this->getInput(scalEq, vectEq);

      // Solve all the linear systems
      this->solveSystems();
      
      // Transfer timestep output back to equations
      this->transferOutput(scalEq, vectEq);
   }

   void SparseLinearCoordinator::addSolverD(const int start)
   {
      SharedSparseDLinearSolver spSolver(new SparseDLinearSolver(start));

      this->mDSolvers.push_back(spSolver);
   }

   void SparseLinearCoordinator::addSolverZ(const int start)
   {
      SharedSparseZLinearSolver spSolver(new SparseZLinearSolver(start));

      this->mZSolvers.push_back(spSolver);
   }

   void SparseLinearCoordinator::buildSolverMatrix(SolverD_iterator solDIt, const int matIdx, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const int idx)
   {
      // Resize matrix if necessary
      if(solDIt->mLHSMatrix.at(matIdx).size() == 0)
      {
         solDIt->mLHSMatrix.at(matIdx).resize(spEq->couplingInfo(comp).systemN(idx), spEq->couplingInfo(comp).systemN(idx));
      }

      solDIt->mLHSMatrix.at(matIdx) += spEq->operatorRow(Equations::IEquation::BOUNDARYROW, comp, idx).first + spEq->operatorRow(Equations::IEquation::LINEARROW, comp, idx).first;
   }

   void SparseLinearCoordinator::buildSolverMatrix(SolverZ_iterator solZIt, const int matIdx, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const int idx, const bool isLhs)
   {
      // Resize matrix if necessary
      if(solZIt->mLHSMatrix.at(matIdx).size() == 0)
      {
         solZIt->mLHSMatrix.at(matIdx).resize(spEq->couplingInfo(comp).systemN(idx), spEq->couplingInfo(comp).systemN(idx));
      }

      DecoupledZSparse bcRow = spEq->operatorRow(Equations::IEquation::BOUNDARYROW, comp, idx);
      DecoupledZSparse linRow = spEq->operatorRow(Equations::IEquation::LINEARROW, comp, idx);
      solZIt->mLHSMatrix.at(matIdx) += linRow.first.cast<MHDComplex>() + MathConstants::cI*linRow.second + bcRow.first.cast<MHDComplex>() + MathConstants::cI*bcRow.second;
   }

}
}
