/** 
 * @file VariableSelector.hpp
 * @brief Definition of some useful typedefs for the variables used in the code
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef VARIABLESELECTOR_HPP
#define VARIABLESELECTOR_HPP

// Configuration includes
//
#include "SmartPointers/SharedPtrMacro.h"

// System includes
//

// External includes
//

// Project includes
//
#include "TypeSelectors/ScalarSelector.hpp"
#include "Variables/Variable.hpp"
#include "Variables/Spectral/ScalarVariable.hpp"
#include "Variables/Spectral/VectorVariable.hpp"

namespace GeoMHDiSCC {

   namespace Datatypes {

      /// Typedef for a ScalarVariable
      typedef Variable<ScalarVariable<SpectralScalarType,2,PhysicalScalarType,3>, 1> ScalarVariableType;

      /// Typedef for a VectorVariable
      typedef Variable<VectorVariable<SpectralScalarType,2,PhysicalScalarType,3>, 1> VectorVariableType;

      /// Typedef for a shared ScalarVariable
      typedef SharedPtrMacro<ScalarVariableType>  SharedScalarVariableType;

      /// Typedef for a shared ScalarVariable
      typedef SharedPtrMacro<VectorVariableType>  SharedVectorVariableType;

   }
}

#endif // VARIABLESELECTOR_HPP
