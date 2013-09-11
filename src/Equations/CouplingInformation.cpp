/** 
 * @file CouplingInformation.cpp
 * @brief Source of the base implementation of a scalar equation
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 * @version 0.9.0
 * @date 2013-09-11
 */

// Debug includes
//
#include "Debug/DebuggerMacro.h"

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
      : mEquationType(TRIVIAL), mHasNonlinear(false), mHasQuasiInverse(false), mHasSource(false), mIsComplex(true), mIndexType(CouplingInformation::SLOWEST), mNSystems(0), mFieldIndex(-1), mSolverIndex(-1), mFieldStart(-1)
   {
   }

   CouplingInformation::~CouplingInformation()
   {
   }

   CouplingInformation::EquationTypeId CouplingInformation::equationType() const
   {
      return this->mEquationType;
   }

   bool CouplingInformation::hasNonlinear() const
   {
      return this->mHasNonlinear;
   }

   bool CouplingInformation::hasQuasiInverse() const
   {
      return this->mHasQuasiInverse;
   }

   bool CouplingInformation::hasSource() const
   {
      return this->mHasSource;
   }

   bool CouplingInformation::isComplex() const
   {
      return this->mIsComplex;
   }

   CouplingInformation::IndexType CouplingInformation::indexType() const
   {
      return this->mIndexType;
   }

   int CouplingInformation::nBlocks() const
   {
      return this->mImplicitFields.size();
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
      return this->blockN(idx)*this->nBlocks();
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

   int CouplingInformation::nExplicit() const
   {
      return this->mExplicitFields.size();
   }

   int CouplingInformation::rhsCols(const int idx) const
   {
      return this->mRhsCols(idx);
   }

   void CouplingInformation::sortImplicitFields(const PhysicalNames::Id fieldId, const FieldComponents::Spectral::Id compId)
   {
      // Sort the implicit fields
      std::sort(this->mImplicitFields.begin(), this->mImplicitFields.end());

      // Extract the position of the equation field
      FieldId_iterator pos = std::find(this->mImplicitFields.begin(), this->mImplicitFields.end(), std::make_pair(fieldId, compId));

      // Set initial field index
      this->mFieldIndex = pos - this->mImplicitFields.begin();

      // Report position of equation field
      DebuggerMacro_showValue("Coupling information field index: ", 1, this->mFieldIndex);
   }

   void CouplingInformation::addImplicitField(const PhysicalNames::Id fieldId, const FieldComponents::Spectral::Id compId)
   {
      this->mImplicitFields.push_back(std::make_pair(fieldId,compId));
   }

   void CouplingInformation::addExplicitField(const PhysicalNames::Id fieldId, const FieldComponents::Spectral::Id compId)
   {
      this->mExplicitFields.push_back(std::make_pair(fieldId,compId));
   }

   void CouplingInformation::setGeneral(const CouplingInformation::EquationTypeId typeId, const bool isComplex, const int fieldStart)
   {
      this->mEquationType = typeId;

      this->mIsComplex = isComplex;

      this->mFieldStart = fieldStart;
   }

   void CouplingInformation::setNonlinear(const bool hasNonlinear, const bool hasQuasiInverse)
   {
      this->mHasNonlinear = hasNonlinear;

      this->mHasQuasiInverse = hasQuasiInverse;
   }

   void CouplingInformation::setSource(const bool hasSource)
   {
      this->mHasSource = hasSource;
   }

   void CouplingInformation::setSizes(const int nSystems, const ArrayI& blockNs, const ArrayI& rhsCols)
   {
      this->mNSystems = nSystems;

      this->mBlockNs = blockNs;

      this->mRhsCols = rhsCols;
   }

   void CouplingInformation::setSolverIndex(const int idx)
   {
      this->mSolverIndex = idx;
   }

   void CouplingInformation::setIndexType(const CouplingInformation::IndexType id)
   {
      this->mIndexType = id;
   }

   CouplingInformation::FieldId_range CouplingInformation::implicitRange() const
   {
      return std::make_pair(this->mImplicitFields.begin(), this->mImplicitFields.end());
   }

   CouplingInformation::FieldId_range CouplingInformation::explicitRange() const
   {
      return std::make_pair(this->mExplicitFields.begin(), this->mExplicitFields.end());
   }
}
}
