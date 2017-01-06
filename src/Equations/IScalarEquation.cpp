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

namespace QuICC {

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

   void IScalarEquation::updateDealiasedUnknown(const Datatypes::SpectralScalarType& rhs, FieldComponents::Spectral::Id compId, Arithmetics::Id arithId)
   {
      // Assert scalar
      assert(compId == FieldComponents::Spectral::SCALAR);

      // Assert dealiasing has taken place!
      assert(this->rUnknown().rDom(0).rPerturbation().data().rows() <= rhs.data().rows());
      assert(this->rUnknown().rDom(0).rPerturbation().data().cols() == rhs.data().cols());

      if(arithId == Arithmetics::SET)
      {
         // Copy values over into unknown
         this->rUnknown().rDom(0).rPerturbation().setData(rhs.data().topRows(this->rUnknown().rDom(0).rPerturbation().data().rows()));

      } else if(arithId == Arithmetics::SETNEG)
      {
         // Copy negative values over into unknown
         this->rUnknown().rDom(0).rPerturbation().setData(-rhs.data().topRows(this->rUnknown().rDom(0).rPerturbation().data().rows()));
      } else if(arithId == Arithmetics::ADD)
      {
         // Add values to unknown
         this->rUnknown().rDom(0).rPerturbation().addData(rhs.data().topRows(this->rUnknown().rDom(0).rPerturbation().data().rows()));
      } else if(arithId == Arithmetics::SUB)
      {
         // Substract values from unknown
         this->rUnknown().rDom(0).rPerturbation().subData(rhs.data().topRows(this->rUnknown().rDom(0).rPerturbation().data().rows()));
      }
   }

   void IScalarEquation::initSpectralMatrices()
   {
      // Make sure it is safe to do nothing
      bool needInit = this->couplingInfo(FieldComponents::Spectral::SCALAR).hasQuasiInverse();

      // Check for Galerkin stencils
      needInit = needInit || this->couplingInfo(FieldComponents::Spectral::SCALAR).isGalerkin();

      // Check explicit linear operators
      CouplingInformation::FieldId_range fRange = this->couplingInfo(FieldComponents::Spectral::SCALAR).explicitRange(ModelOperator::EXPLICIT_LINEAR);
      needInit = needInit || (std::distance(fRange.first,fRange.second) > 0);
      // Check explicit nonlinear operators
      fRange = this->couplingInfo(FieldComponents::Spectral::SCALAR).explicitRange(ModelOperator::EXPLICIT_NONLINEAR);
      needInit = needInit || (std::distance(fRange.first,fRange.second) > 0);
      // Check explicit nextstep operators
      fRange = this->couplingInfo(FieldComponents::Spectral::SCALAR).explicitRange(ModelOperator::EXPLICIT_NEXTSTEP);
      needInit = needInit || (std::distance(fRange.first,fRange.second) > 0);

      // Initialise spectral matrices
      if(needInit)
      {
         this->initSpectralMatricesComponent(this->mspBcIds, FieldComponents::Spectral::SCALAR);
      }
   }

   void IScalarEquation::defineCoupling(FieldComponents::Spectral::Id comp, CouplingInformation::EquationTypeId eqType, const int iZero, const bool hasNL, const bool hasSource, const bool allowExplicit)
   {
      this->dispatchCoupling(comp, eqType, iZero, hasNL, hasSource, this->unknown().dom(0).spRes(), allowExplicit);
   }

   void  IScalarEquation::buildModelMatrix(DecoupledZSparse& rModelMatrix, const ModelOperator::Id opId, FieldComponents::Spectral::Id comp, const int matIdx, const ModelOperatorBoundary::Id bcType) const
   {
      this->dispatchModelMatrix(rModelMatrix, opId, comp, matIdx, bcType, this->unknown().dom(0).spRes(), this->couplingInfo(FieldComponents::Spectral::SCALAR).eigenTools().getEigs(this->spRes(), matIdx));
   }

   void IScalarEquation::setGalerkinStencil(FieldComponents::Spectral::Id comp, SparseMatrix &mat, const int matIdx) const
   {
      this->dispatchGalerkinStencil(comp, mat, matIdx, this->unknown().dom(0).spRes(), this->couplingInfo(FieldComponents::Spectral::SCALAR).eigenTools().getEigs(this->spRes(), matIdx));
   }

   void IScalarEquation::setExplicitBlock(FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const ModelOperator::Id opId, const SpectralFieldId fieldId, const int matIdx) const
   {
      this->dispatchExplicitBlock(compId, mat, opId, fieldId, matIdx, this->unknown().dom(0).spRes(), this->couplingInfo(FieldComponents::Spectral::SCALAR).eigenTools().getEigs(this->spRes(), matIdx));
   }

   void IScalarEquation::setNLComponents()
   {
      this->addNLComponent(FieldComponents::Spectral::SCALAR, 0);
   }
}
}
