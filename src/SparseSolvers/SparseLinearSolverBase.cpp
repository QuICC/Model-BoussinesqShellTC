/** \file SparseLinearSolverBase.cpp
 *  \brief Implementation of the base for linear solver structures
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "SparseSolvers/SparseLinearSolverBase.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Solver {

   SparseLinearSolverBase::SparseLinearSolverBase(const int start)
      : mZeroIdx(start)
   {
      // Safety assert
      assert(start >= 0);
   }

   SparseLinearSolverBase::~SparseLinearSolverBase()
   {
   }

   void SparseLinearSolverBase::addInformation(const SpectralFieldId& id, const ArrayI& startRow)
   {
      this->mFieldIds.push_back(id);

      this->mInformation.insert(std::make_pair(id, startRow));
   }

   int SparseLinearSolverBase::startRow(const SpectralFieldId& id, const int i) const
   {
      return this->mInformation.find(id)->second(i);
   }

   SparseLinearSolverBase::FieldId_range SparseLinearSolverBase::fieldRange() const
   {
      return std::make_pair(this->mFieldIds.begin(), this->mFieldIds.end());
   }
}
}
