/** \file CartesianSpectralOperatorTypedefs.hpp
 *  \brief Definition of some useful typedefs for the variables used in the cartesian code
 */

#ifndef CARTESIANSPECTRALOPERATORTYPEDEFS_HPP
#define CARTESIANSPECTRALOPERATORTYPEDEFS_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "Base/SpectralOperators/UnitSpectralOperator.hpp"
#include "Base/SpectralOperators/ChebyshevSpectralOperator.hpp"

namespace EPMPhoenix {

   namespace Code {

      #if defined EPMPHOENIX_CODEDIM_1D
         #error "Cartesian code in 1D is not implemented"

      #elif defined EPMPHOENIX_CODEDIM_2D
         #if defined EPMPHOENIX_SPATIALSCHEME_TF
            /// Typedef for the first dimension
            typedef ChebyshevSpectralOperator SpectralOperator1DType;
         #elif defined EPMPHOENIX_SPATIALSCHEME_TPF
            /// Typedef for the first dimension
            typedef ChebyshevSpectralOperator SpectralOperator1DType;
         #endif // EPMPHOENIX_SPATIALSCHEME_TF

         /// Typedef for the second dimension
         typedef UnitSpectralOperator SpectralOperator2DType;

      #elif defined EPMPHOENIX_CODEDIM_3D
         #if defined EPMPHOENIX_SPATIALSCHEME_TFT
            /// Typedef for the first dimension
            typedef ChebyshevSpectralOperator SpectralOperator1DType;

            /// Typedef for the third dimension
            typedef ChebyshevSpectralOperator SpectralOperator3DType;

         #elif defined EPMPHOENIX_SPATIALSCHEME_TPFTP
            /// Typedef for the first dimension
            typedef ChebyshevSpectralOperator SpectralOperator1DType;

            /// Typedef for the third dimension
            typedef ChebyshevSpectralOperator SpectralOperator3DType;

         #elif defined EPMPHOENIX_SPATIALSCHEME_TFF
            /// Typedef for the first dimension
            typedef ChebyshevSpectralOperator SpectralOperator1DType;

            /// Typedef for the third dimension
            typedef UnitSpectralOperator SpectralOperator3DType;

         #elif defined EPMPHOENIX_SPATIALSCHEME_TPFF
            /// Typedef for the first dimension
            typedef ChebyshevSpectralOperator SpectralOperator1DType;

            /// Typedef for the third dimension
            typedef UnitSpectralOperator SpectralOperator3DType;

         #endif // EPMPHOENIX_SPATIALSCHEME_TF

         /// Typedef for the second dimension
         typedef UnitSpectralOperator SpectralOperator2DType;

      #else
         #error "Number of dimension of simulation code not specified"
      #endif // EPMPHOENIX_CODEDIM_1D

   }
}

#endif // CARTESIANSPECTRALOPERATORTYPEDEFS_HPP
