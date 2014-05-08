/** 
 * @file SparseSolverBase.cpp
 * @brief Implementation of the base for linear solver structures
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
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
#include "Exceptions/Exception.hpp"

namespace GeoMHDiSCC {

namespace Solver {

   SparseSolverBase::SparseSolverBase(const int start)
      : mZeroIdx(start), mPyName(""), mIsInitialized(false)
   {
      // Safety assert
      assert(start >= 0);
   }

   SparseSolverBase::~SparseSolverBase()
   {
   }

   void SparseSolverBase::addInformation(const SpectralFieldId& id, const ArrayI& startRow, const std::string& pyName)
   {
      if(this->mPyName.length() == 0)
      {
         this->mPyName = pyName;
      } else if(this->mPyName.compare(pyName) != 0)
      {
         throw Exception("Attempted to setup solver with different python generators!");
      }

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

   bool SparseSolverBase::isInitialized() const
   {
      return this->mIsInitialized;
   }

   void SparseSolverBase::setInitialized()
   {
      this->mIsInitialized = true;
   }
}
}
