/** \file SparseLinearSolverBaseBase.cpp
 *  \brief Implementation of the base for linear solver structures
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "SparseSolvers/SparseLinearSolverBaseBase.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Solver {

   SparseLinearSolverBaseBase::SparseLinearSolverBaseBase(const int start)
      : mZeroIdx(start)
   {
      // Safety assert
      assert(start >= 0);
   }

   SparseLinearSolverBaseBase::~SparseLinearSolverBaseBase()
   {
   }

   void SparseLinearSolverBaseBase::addInformation(const SpectralFieldId& id, const ArrayI& startRow)
   {
      this->mFieldIds.push_back(id);

      this->mInformation.insert(std::make_pair(id, startRow));
   }

   int SparseLinearSolverBaseBase::startRow(const SpectralFieldId& id, const int i) const
   {
      return this->mInformation.find(id)->second(i);
   }

   SparseLinearSolverBaseBase::FieldId_range SparseLinearSolverBaseBase::fieldRange() const
   {
      return std::make_pair(this->mFieldIds.begin(), this->mFieldIds.end());
   }
}
}
