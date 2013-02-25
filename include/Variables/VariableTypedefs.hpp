/** \file CartesianVariableTypedefs.hpp
 *  \brief Definition of some useful typedefs for the variables used in the cartesian code
 */

#ifndef CARTESIANVARIABLETYPEDEFS_HPP
#define CARTESIANVARIABLETYPEDEFS_HPP

// Configuration includes
//
#include "Base/PrepMacros/SmartPointerMacro.h"

// System includes
//

// External includes
//

// Project includes
//
#include "Cartesian/General/CartesianScalarTypedefs.hpp"
#include "Simulation/Variables/Variable.hpp"
#include "Simulation/Variables/Spectral/ScalarVariable.hpp"
#include "Simulation/Variables/Spectral/VectorVariable.hpp"

namespace EPMPhoenix {

   namespace Code {

      #if defined EPMPHOENIX_CODEDIM_1D
         #error "Cartesian code in 1D is not possible"
      #elif defined EPMPHOENIX_CODEDIM_2D
         /// Typedef for a ScalarVariable
         typedef EPMPhoenix::Variable<EPMPhoenix::ScalarVariable<SpectralScalarType,2,PhysicalScalarType,2>, 1> ScalarVariable;

         /// Typedef for a VectorVariable
         typedef EPMPhoenix::Variable<EPMPhoenix::VectorVariable<SpectralScalarType,2,PhysicalScalarType,2>, 1> VectorVariable;
      #elif defined EPMPHOENIX_CODEDIM_3D
         /// Typedef for a ScalarVariable
         typedef EPMPhoenix::Variable<EPMPhoenix::ScalarVariable<SpectralScalarType,2,PhysicalScalarType,3>, 1> ScalarVariable;

         /// Typedef for a VectorVariable
         typedef EPMPhoenix::Variable<EPMPhoenix::VectorVariable<SpectralScalarType,2,PhysicalScalarType,3>, 1> VectorVariable;
      #else
         #error "Number of dimension is simulation code not specified"
      #endif // EPMPHOENIX_CODEDIM_1D

      /// Typedef for a shared ScalarVariable
      typedef SharedPtrMacro<ScalarVariable>  SharedScalarVariable;

      /// Typedef for a shared ScalarVariable
      typedef SharedPtrMacro<VectorVariable>  SharedVectorVariable;

   }
}

#endif // CARTESIANVARIABLETYPEDEFS_HPP
