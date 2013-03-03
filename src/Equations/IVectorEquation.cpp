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

   IVectorEquation::IVectorEquation(SharedIEquationParameters spEqParams)
      : IEvolutionEquation(spEqParams)
   {
   }

   IVectorEquation::~IVectorEquation()
   {
   }

   void IVectorEquation::setUnknown(Datatypes::SharedVectorVariableType spUnknown)
   {
      /// \mhdBug Fake implementation
      
      this->mspUnknown = spUnknown;
   }

   const Datatypes::VectorVariableType& IVectorEquation::unknown() const
   {
      /// \mhdBug Fake implementation
      
      return *this->mspUnknown;
   }

   Datatypes::VectorVariableType& IVectorEquation::rUnknown()
   {
      /// \mhdBug Fake implementation
      
      return *this->mspUnknown;
   }

   void IVectorEquation::computeLinear(Datatypes::SpectralScalarType& rRHS, FieldComponents::Spectral::Id id)
   {
      /// \mhdBug Fake implementation
   }
}
