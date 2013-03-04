/** \file IVectorEquation.cpp
 *  \brief Source of the base implementation of a vector equation
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

namespace GeoMHDiSCC {

namespace Equations {

   IVectorEquation::IVectorEquation(SharedIEquationParameters spEqParams)
      : IEvolutionEquation(spEqParams)
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

   void IVectorEquation::computeLinear(Datatypes::SpectralScalarType& rRHS, FieldComponents::Spectral::Id id) const
   {
      // Empty default implementation
   }

   void IVectorEquation::prepareTimestep(const Datatypes::SpectralScalarType& rhs, FieldComponents::Spectral::Id id)
   {
      // Assert dealiasing has taken place!
      assert(this->rUnknown().rDom(0).rPerturbation().rComp(id).data().rows() < rhs.data().rows());
      assert(this->rUnknown().rDom(0).rPerturbation().rComp(id).data().cols() < rhs.data().cols());

      // Copy values over into unknown
      this->rUnknown().rDom(0).rPerturbation().rComp(id).setData(rhs.data().topRows(this->rUnknown().rDom(0).rPerturbation().rComp(id).data().rows()));
   }
}
}
