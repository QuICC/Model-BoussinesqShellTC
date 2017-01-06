/** 
 * @file EquationConditions.cpp
 * @brief Source of the equation condition functors
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// System includes
//
#include <cassert>
#include <algorithm>

// External includes
//

// Class include
//
#include "Equations/Tools/EquationConditions.hpp"

// Project includes
//

namespace QuICC {

namespace Equations {

namespace Conditions {

   bool IsPrognostic::operator()(SharedIScalarEquation eqA)
   {
      return eqA->couplingInfo(FieldComponents::Spectral::SCALAR).equationType() == CouplingInformation::PROGNOSTIC;
   }

   bool IsPrognostic::operator()(SharedIVectorEquation eqA)
   {
      return eqA->couplingInfo(FieldComponents::Spectral::ONE).equationType() == CouplingInformation::PROGNOSTIC;
   }

   bool IsDiagnostic::operator()(SharedIScalarEquation eqA)
   {
      return eqA->couplingInfo(FieldComponents::Spectral::SCALAR).equationType() == CouplingInformation::DIAGNOSTIC;
   }

   bool IsDiagnostic::operator()(SharedIVectorEquation eqA)
   {
      return eqA->couplingInfo(FieldComponents::Spectral::ONE).equationType() == CouplingInformation::DIAGNOSTIC;
   }

   bool IsTrivial::operator()(SharedIScalarEquation eqA)
   {
      return eqA->couplingInfo(FieldComponents::Spectral::SCALAR).equationType() == CouplingInformation::TRIVIAL;
   }

   bool IsTrivial::operator()(SharedIVectorEquation eqA)
   {
      return eqA->couplingInfo(FieldComponents::Spectral::ONE).equationType() == CouplingInformation::TRIVIAL;
   }

   bool IsWrapper::operator()(SharedIScalarEquation eqA)
   {
      return eqA->couplingInfo(FieldComponents::Spectral::SCALAR).equationType() == CouplingInformation::WRAPPER;
   }

   bool IsWrapper::operator()(SharedIVectorEquation eqA)
   {
      return eqA->couplingInfo(FieldComponents::Spectral::ONE).equationType() == CouplingInformation::WRAPPER;
   }
}
}
}
