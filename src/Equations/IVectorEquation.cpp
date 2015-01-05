/** 
 * @file IVectorEquation.cpp
 * @brief Source of the base implementation of a vector equation
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
#include "Equations/IVectorEquation.hpp"

// Project includes
//
#include "TypeSelectors/EquationEigenSelector.hpp"
#include "Exceptions/Exception.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   IVectorEquation::IVectorEquation(SharedEquationParameters spEqParams)
      : IEquation(spEqParams)
   {
   }

   IVectorEquation::~IVectorEquation()
   {
   }

   void IVectorEquation::setUnknown(Datatypes::SharedVectorVariableType spUnknown)
   {
      this->mspUnknown = spUnknown;
   }

   const Datatypes::VectorVariableType& IVectorEquation::unknown() const
   {
      // Safety assert
      assert(this->mspUnknown);

      return *this->mspUnknown;
   }

   Datatypes::VectorVariableType& IVectorEquation::rUnknown()
   {
      // Safety assert
      assert(this->mspUnknown);

      return *this->mspUnknown;
   }

   SharedResolution IVectorEquation::spRes() const
   {
      return this->unknown().dom(0).spRes();
   }

   int IVectorEquation::nSpectral() const
   {
      return this->mRequirements.field(this->name()).spectralIds().size();
   }

   IVectorEquation::SpectralComponent_range IVectorEquation::spectralRange() const
   {
      return std::make_pair(this->mRequirements.field(this->name()).spectralIds().begin(), this->mRequirements.field(this->name()).spectralIds().end());
   }

   void IVectorEquation::updateDealiasedUnknown(const Datatypes::SpectralScalarType& rhs, FieldComponents::Spectral::Id compId, Arithmetics::Id arithId)
   {
      // Assert dealiasing has taken place!
      assert(this->rUnknown().rDom(0).rPerturbation().rComp(compId).data().rows() < rhs.data().rows());
      assert(this->rUnknown().rDom(0).rPerturbation().rComp(compId).data().cols() == rhs.data().cols());

      if(arithId == Arithmetics::SET)
      {
         // Copy values over into unknown
         this->rUnknown().rDom(0).rPerturbation().rComp(compId).setData(rhs.data().topRows(this->rUnknown().rDom(0).rPerturbation().rComp(compId).data().rows()));
      } else if(arithId == Arithmetics::SETNEG)
      {
         // Copy negative values over into unknown
         this->rUnknown().rDom(0).rPerturbation().rComp(compId).setData(-rhs.data().topRows(this->rUnknown().rDom(0).rPerturbation().rComp(compId).data().rows()));
      } else if(arithId == Arithmetics::ADD)
      {
         // Add values to unknown
         this->rUnknown().rDom(0).rPerturbation().rComp(compId).addData(rhs.data().topRows(this->rUnknown().rDom(0).rPerturbation().rComp(compId).data().rows()));
      } else if(arithId == Arithmetics::SUB)
      {
         // Substract values from unknown
         this->rUnknown().rDom(0).rPerturbation().rComp(compId).subData(rhs.data().topRows(this->rUnknown().rDom(0).rPerturbation().rComp(compId).data().rows()));
      }
   }

   void IVectorEquation::initSpectralMatrices(const SharedSimulationBoundary spBcIds)
   {
      // Store the boundary condition list
      this->mspBcIds = spBcIds;

      IVectorEquation::SpectralComponent_range range = this->spectralRange();

      for(SpectralComponent_iterator it = range.first; it != range.second; ++it)
      {
         // Make sure it is safe to do nothing
         bool needInit = this->couplingInfo(*it).hasQuasiInverse();

         CouplingInformation::FieldId_range fRange = this->couplingInfo(*it).explicitRange();
         needInit = needInit || (std::distance(fRange.first, fRange.second) > 0);

         // Initialise spectral matrices
         if(needInit)
         {
            this->initSpectralMatricesComponent(spBcIds, *it);
         }
      }
   }

   void IVectorEquation::defineCoupling(FieldComponents::Spectral::Id comp, CouplingInformation::EquationTypeId eqType, const int iZero, const bool hasNL, const bool hasQI, const bool hasSource, const bool allowExplicit)
   {
      this->dispatchCoupling(comp, eqType, iZero, hasNL, hasQI, hasSource, this->unknown().dom(0).spRes(), allowExplicit);
   }

   void  IVectorEquation::buildModelMatrix(DecoupledZSparse& rModelMatrix, const ModelOperator::Id opId, FieldComponents::Spectral::Id comp, const int matIdx, const ModelOperatorBoundary::Id bcType) const
   {
      this->dispatchModelMatrix(rModelMatrix, opId, comp, matIdx, bcType, this->unknown().dom(0).spRes(), EigenSelector::getEigs(*this, matIdx));
   }

   void IVectorEquation::setGalerkinStencil(FieldComponents::Spectral::Id comp, SparseMatrix &mat, const int matIdx) const
   {
      this->dispatchGalerkinStencil(comp, mat, matIdx, this->unknown().dom(0).spRes(), EigenSelector::getEigs(*this, matIdx));
   }

   void IVectorEquation::setQuasiInverse(FieldComponents::Spectral::Id comp, SparseMatrix &mat, const int matIdx) const
   {
      this->dispatchQuasiInverse(comp, mat, matIdx, this->unknown().dom(0).spRes(), EigenSelector::getEigs(*this, matIdx));
   }

   void IVectorEquation::setExplicitLinearBlock(FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const int matIdx) const
   {
      this->dispatchExplicitLinearBlock(compId, mat, fieldId, matIdx, this->unknown().dom(0).spRes(), EigenSelector::getEigs(*this, matIdx));
   }
}
}
