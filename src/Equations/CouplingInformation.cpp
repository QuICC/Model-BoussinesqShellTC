/** 
 * @file CouplingInformation.cpp
 * @brief Source of the base implementation of a scalar equation
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
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
#include "TypeSelectors/EquationEigenSelector.hpp"

namespace QuICC {

namespace Equations {

   CouplingInformation::CouplingInformation()
      : mEquationType(TRIVIAL), mHasNonlinear(false), mHasQuasiInverse(false), mHasSource(false), mIsComplex(true), mIsGalerkin(false), mIndexType(CouplingInformation::SLOWEST_SINGLE_RHS), mNSystems(0), mFieldIndex(-1), mSolverIndex(-1), mFieldStart(-1)
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

   bool CouplingInformation::isGalerkin() const
   {
      return this->mIsGalerkin;
   }

   CouplingInformation::IndexType CouplingInformation::indexType() const
   {
      return this->mIndexType;
   }

   const IEigenTools& CouplingInformation::eigenTools() const
   {
      return *(this->mspEigenTools);
   }

   int CouplingInformation::nBlocks() const
   {
      return this->mImplicitFields.size();
   }

   int CouplingInformation::nSystems() const
   {
      return this->mNSystems;
   }

   int CouplingInformation::tauN(const int idx) const
   {
      return this->mTauNs(idx);
   }

   int CouplingInformation::galerkinN(const int idx) const
   {
      return this->mGalerkinNs(idx);
   }

   int CouplingInformation::galerkinShift(const int idx, const int dim) const
   {
      return this->mGalerkinShifts(idx, dim);
   }

   int CouplingInformation::systemN(const int idx) const
   {
      return this->mSystemNs(idx);
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

   void CouplingInformation::sortImplicitFields(const PhysicalNames::Id fieldId, const FieldComponents::Spectral::Id compId)
   {
      // Sort the implicit fields
      std::sort(this->mImplicitFields.begin(), this->mImplicitFields.end());

      // Extract the position of the equation field
      FieldId_iterator pos = std::find(this->mImplicitFields.begin(), this->mImplicitFields.end(), std::make_pair(fieldId, compId));
      assert((this->equationType() == TRIVIAL && pos == this->mImplicitFields.end()) || (this->equationType() == WRAPPER && pos == this->mImplicitFields.end()) || pos != this->mImplicitFields.end());

      // Set initial field index
      this->mFieldIndex = pos - this->mImplicitFields.begin();

      // Report position of equation field
      DebuggerMacro_showValue("Coupling information field index: ", 1, this->mFieldIndex);
   }

   void CouplingInformation::addImplicitField(const PhysicalNames::Id fieldId, const FieldComponents::Spectral::Id compId)
   {
      this->mImplicitFields.push_back(std::make_pair(fieldId,compId));
   }

   void CouplingInformation::addExplicitField(const PhysicalNames::Id fieldId, const FieldComponents::Spectral::Id compId, const ModelOperator::Id opId)
   {
      if(opId == ModelOperator::EXPLICIT_LINEAR)
      {
         this->mExplicitLFields.push_back(std::make_pair(fieldId,compId));

      } else if(opId == ModelOperator::EXPLICIT_NONLINEAR)
      {
         this->mExplicitNLFields.push_back(std::make_pair(fieldId,compId));

      } else if(opId == ModelOperator::EXPLICIT_NEXTSTEP)
      {
         this->mExplicitNSFields.push_back(std::make_pair(fieldId,compId));
      }

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

   void CouplingInformation::setSizes(const int nSystems, const ArrayI& tauNs, const ArrayI& galerkinNs, const MatrixI& galerkinShifts, const ArrayI& rhsCols, const ArrayI& systemNs)
   {
      this->mNSystems = nSystems;

      this->mTauNs = tauNs;

      this->mGalerkinNs = galerkinNs;

      this->mGalerkinShifts = galerkinShifts;

      this->mRhsCols = rhsCols;

      this->mIsGalerkin = (this->mGalerkinShifts.sum() > 0);

      this->mSystemNs = systemNs;
   }

   void CouplingInformation::setSolverIndex(const int idx)
   {
      this->mSolverIndex = idx;
   }

   void CouplingInformation::setIndexType(const CouplingInformation::IndexType id)
   {
      this->mIndexType = id;

      this->mspEigenTools = eigenSelector(this->mIndexType);
   }

   CouplingInformation::FieldId_range CouplingInformation::implicitRange() const
   {
      return std::make_pair(this->mImplicitFields.begin(), this->mImplicitFields.end());
   }

   CouplingInformation::FieldId_range CouplingInformation::explicitRange(const ModelOperator::Id opId) const
   {
      CouplingInformation::FieldId_range range;
      if(opId == ModelOperator::EXPLICIT_LINEAR)
      {
         range = std::make_pair(this->mExplicitLFields.begin(), this->mExplicitLFields.end());

      } else if(opId == ModelOperator::EXPLICIT_NONLINEAR)
      {
         range = std::make_pair(this->mExplicitNLFields.begin(), this->mExplicitNLFields.end());

      } else if(opId == ModelOperator::EXPLICIT_NEXTSTEP)
      {
         range = std::make_pair(this->mExplicitNSFields.begin(), this->mExplicitNSFields.end());
      }

      return range;
   }
}
}
