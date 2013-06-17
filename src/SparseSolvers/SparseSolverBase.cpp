/** \file SparseSolverBase.cpp
 *  \brief Implementation of the base for linear solver structures
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "SparseSolvers/SparseSolverBase.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Solver {

   SparseSolverBase::SparseSolverBase(const int start)
      : mZeroIdx(start)
   {
      // Safety assert
      assert(start >= 0);
   }

   SparseSolverBase::~SparseSolverBase()
   {
   }

   void SparseSolverBase::addInformation(const SpectralFieldId& id, const ArrayI& startRow)
   {
      this->mFieldIds.push_back(id);

      this->mInformation.insert(std::make_pair(id, startRow));
   }

   int SparseSolverBase::startRow(const SpectralFieldId& id, const int i) const
   {
      return this->mInformation.find(id)->second(i);
   }

   SparseSolverBase::FieldId_range SparseSolverBase::fieldRange() const
   {
      return std::make_pair(this->mFieldIds.begin(), this->mFieldIds.end());
   }
}
}
