/** \file IVectorPEquation.cpp
 *  \brief Source of the base implementation of a vector prognostic equation
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "Equations/IVectorPEquation.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Equations {

   IVectorPEquation::IVectorPEquation(SharedEquationParameters spEqParams)
      : IPrognosticEquation(spEqParams)
   {
   }

   IVectorPEquation::~IVectorPEquation()
   {
   }

   void IVectorPEquation::setUnknown(Datatypes::SharedVectorVariableType spUnknown)
   {
      this->mspUnknown = spUnknown;
   }

   const Datatypes::VectorVariableType& IVectorPEquation::unknown() const
   {
      // Safety assert
      assert(this->mspUnknown);

      return *this->mspUnknown;
   }

   Datatypes::VectorVariableType& IVectorPEquation::rUnknown()
   {
      // Safety assert
      assert(this->mspUnknown);

      return *this->mspUnknown;
   }

   void IVectorPEquation::prepareTimestep(const Datatypes::SpectralScalarType& rhs, FieldComponents::Spectral::Id id)
   {
      // Assert dealiasing has taken place!
      assert(this->rUnknown().rDom(0).rPerturbation().rComp(id).data().rows() < rhs.data().rows());
      assert(this->rUnknown().rDom(0).rPerturbation().rComp(id).data().cols() < rhs.data().cols());

      // Copy values over into unknown
      this->rUnknown().rDom(0).rPerturbation().rComp(id).setData(rhs.data().topRows(this->rUnknown().rDom(0).rPerturbation().rComp(id).data().rows()));
   }
}
}
