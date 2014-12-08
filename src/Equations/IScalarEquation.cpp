/** 
 * @file IScalarEquation.cpp
 * @brief Source of the base implementation of a scalar equation
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "Equations/IScalarEquation.hpp"

// Project includes
//
#include "Base/MathConstants.hpp"
#include "TypeSelectors/EquationEigenSelector.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   IScalarEquation::IScalarEquation(SharedEquationParameters spEqParams)
      : IEquation(spEqParams)
   {
   }

   IScalarEquation::~IScalarEquation()
   {
   }

   void IScalarEquation::setUnknown(Datatypes::SharedScalarVariableType spUnknown)
   {
      this->mspUnknown = spUnknown;
   }

   const Datatypes::ScalarVariableType& IScalarEquation::unknown() const
   {
      // Safety assert
      assert(this->mspUnknown);

      return *this->mspUnknown;
   }

   Datatypes::ScalarVariableType& IScalarEquation::rUnknown()
   {
      // Safety assert
      assert(this->mspUnknown);

      return *this->mspUnknown;
   }

   SharedResolution IScalarEquation::spRes() const
   {
      return this->unknown().dom(0).spRes();
   }

   void IScalarEquation::updateDealiasedUnknown(const Datatypes::SpectralScalarType& rhs, FieldComponents::Spectral::Id compId)
   {
      // Assert scalar
      assert(compId == FieldComponents::Spectral::SCALAR);

      // Assert dealiasing has taken place!
      assert(this->rUnknown().rDom(0).rPerturbation().data().rows() < rhs.data().rows());
      assert(this->rUnknown().rDom(0).rPerturbation().data().cols() == rhs.data().cols());

      // Copy values over into unknown
      this->rUnknown().rDom(0).rPerturbation().setData(rhs.data().topRows(this->rUnknown().rDom(0).rPerturbation().data().rows()));
   }

   void IScalarEquation::initSpectralMatrices(const SharedSimulationBoundary spBcIds)
   {
      // Store the boundary condition list
      this->mspBcIds = spBcIds;

      // Make sure it is safe to do nothing
      bool needInit = this->couplingInfo(FieldComponents::Spectral::SCALAR).hasQuasiInverse();

      CouplingInformation::FieldId_range fRange = this->couplingInfo(FieldComponents::Spectral::SCALAR).explicitRange();
      needInit = needInit || (std::distance(fRange.first,fRange.second) > 0);

      // Initialise spectral matrices
      if(needInit)
      {
         this->initSpectralMatricesComponent(spBcIds, FieldComponents::Spectral::SCALAR);
      }
   }

   void IScalarEquation::defineCoupling(FieldComponents::Spectral::Id comp, CouplingInformation::EquationTypeId eqType, const int iZero, const bool hasNL, const bool hasQI, const bool hasSource, const bool allowExplicit)
   {
      this->dispatchCoupling(comp, eqType, iZero, hasNL, hasQI, hasSource, this->unknown().dom(0).spRes(), allowExplicit);
   }

   void  IScalarEquation::buildModelMatrix(DecoupledZSparse& rModelMatrix, const ModelOperator::Id opId, FieldComponents::Spectral::Id comp, const int matIdx, const ModelOperatorBoundary::Id bcType) const
   {
      this->dispatchModelMatrix(rModelMatrix, opId, comp, matIdx, bcType, this->unknown().dom(0).spRes(), EigenSelector::getEigs(*this, matIdx));
   }

   void IScalarEquation::setGalerkinStencil(FieldComponents::Spectral::Id comp, SparseMatrix &mat, const int matIdx) const
   {
      this->dispatchGalerkinStencil(comp, mat, matIdx, this->unknown().dom(0).spRes(), EigenSelector::getEigs(*this, matIdx));
   }

   void IScalarEquation::setQuasiInverse(FieldComponents::Spectral::Id comp, SparseMatrix &mat, const int matIdx) const
   {
      this->dispatchQuasiInverse(comp, mat, matIdx, this->unknown().dom(0).spRes(), EigenSelector::getEigs(*this, matIdx));
   }

   void IScalarEquation::setExplicitLinearBlock(FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const int matIdx) const
   {
      this->dispatchExplicitLinearBlock(compId, mat, fieldId, matIdx, this->unknown().dom(0).spRes(), EigenSelector::getEigs(*this, matIdx));
   }
}
}
