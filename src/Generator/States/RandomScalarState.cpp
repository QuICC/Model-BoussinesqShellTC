/** \file RandomScalarState.cpp
 *  \brief Source of the implementation of the general random scalar state equation
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "Generator/States/RandomScalarState.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "TypeSelectors/SpectralSelector.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   RandomScalarState::RandomScalarState(SharedEquationParameters spEqParams, const PhysicalNames::Id name)
      : IScalarEquation(spEqParams)
   {
      // Set name of unknown
      this->setName(name);

      // Set the variable requirements
      this->setRequirements();
   }

   RandomScalarState::~RandomScalarState()
   {
   }

   void RandomScalarState::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const
   {
      // Assert on scalar component is used
      assert(id == FieldComponents::Physical::SCALAR);

   }

   DecoupledZSparse RandomScalarState::operatorRow(const IEquation::OperatorRowId opId, FieldComponents::Spectral::Id compId, const int matIdx) const
   {
      if(opId == IEquation::LINEARROW)
      {
         return linearRow1DPeriodic(*this, compId, matIdx);
      } else if(opId == IEquation::BOUNDARYROW)
      {
         return boundaryRow1DPeriodic(*this, compId, matIdx);
      } else
      {
         throw Exception("Unknown operator row ID");
      }
   }

   void RandomScalarState::setRequirements()
   {
      // Add unknown to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(this->name(), FieldRequirement(true, true, false, false));
   }

   void RandomScalarState::setCoupling()
   {
   }

   void RandomScalarState::setQuasiInverse(SparseMatrix& mat) const
   {
      Equations::quasiInverse(*this, mat);
   }

   void RandomScalarState::setExplicitLinearBlock(DecoupledZSparse& mat, const SpectralFieldId fieldId, const MHDFloat k) const
   {
      Equations::linearBlock(*this, mat, fieldId, k);
   }

}
}
