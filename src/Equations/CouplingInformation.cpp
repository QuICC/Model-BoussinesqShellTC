/** \file CouplingInformation.cpp
 *  \brief Source of the base implementation of a scalar equation
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "Equations/CouplingInformation.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Equations {

   CouplingInformation::CouplingInformation()
      : mIsComplex(true), mNSystems(0), mFieldIndex(-1), mSolverIndex(-1), mFieldStart(-1)
   {
   }

   CouplingInformation::~CouplingInformation()
   {
   }

   bool CouplingInformation::isComplex() const
   {
      return this->mIsComplex;
   }

   int CouplingInformation::nBlocks() const
   {
      return this->mFields.size();
   }

   int CouplingInformation::nSystems() const
   {
      return this->mNSystems;
   }

   int CouplingInformation::blockN(const int idx) const
   {
      return this->mBlockNs(idx);
   }

   int CouplingInformation::systemN(const int idx) const
   {
      return this->blockN(idx)*this->mFields.size();
   }

   int CouplingInformation::fieldIndex() const
   {
      return this->mFieldIndex;
   }

   int CouplingInformation::solverIndex() const
   {
      return this->mSolverIndex;
   }

   int CouplingInformation::fieldStart() const
   {
      return this->mFieldStart;
   }

   int CouplingInformation::rhsCols(const int idx) const
   {
      return this->mRhsCols(idx);
   }

   void CouplingInformation::addField(const PhysicalNames::Id fieldId, const FieldComponents::Spectral::Id compId, const bool isSelf)
   {
      // If field is itself, set the field index
      if(isSelf)
      {
         this->mFieldIndex = this->mFields.size();
      }

      this->mFields.push_back(std::make_pair(fieldId,compId));
   }

   void CouplingInformation::setGeneral(const int solverIndex, const bool isComplex, const int fieldStart)
   {
      this->mSolverIndex = solverIndex;

      this->mIsComplex = isComplex;

      this->mFieldStart = fieldStart;
   }

   void CouplingInformation::setSizes(const int nSystems, const ArrayI& blockNs, const ArrayI& rhsCols)
   {
      this->mNSystems = nSystems;

      this->mBlockNs = blockNs;

      this->mRhsCols = rhsCols;
   }

   CouplingInformation::field_iterator_range CouplingInformation::fieldRange() const
   {
      return std::make_pair(this->mFields.begin(), this->mFields.end());
   }
}
}
