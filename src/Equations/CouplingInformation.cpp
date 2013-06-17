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

   void CouplingInformation::addImplicitField(const PhysicalNames::Id fieldId, const FieldComponents::Spectral::Id compId, const bool isSelf)
   {
      // If field is itself, set the field index
      if(isSelf)
      {
         this->mFieldIndex = this->mImplicitFields.size();
      }

      this->mImplicitFields.push_back(std::make_pair(fieldId,compId));
   }

   void CouplingInformation::addExplicitField(const PhysicalNames::Id fieldId, const FieldComponents::Spectral::Id compId)
   {
      this->mExplicitFields.push_back(std::make_pair(fieldId,compId));
   }

   void CouplingInformation::setGeneral(const int solverIndex, const bool isComplex, const int fieldStart)
   {
      // Timestepping is required
      if(solverIndex > 0)
      {
         this->mEquationType = PROGNOSTIC;
         this->mSolverIndex = solverIndex - 1;

      // Solver is required (but no time marching)
      } else if(solverIndex < 0)
      {
         this->mEquationType = DIAGNOSTIC;
         this->mSolverIndex = std::abs(solverIndex) - 1;

      // No solver and no timestepping
      } else
      {
         this->mEquationType = TRIVIAL;
         this->mSolverIndex = 0;
      }

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
