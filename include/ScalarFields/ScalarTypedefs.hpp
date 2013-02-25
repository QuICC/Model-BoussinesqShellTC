/** \file ScalarTypedefs.hpp
 *  \brief Definition of some useful typedefs for the scalars used in the cartesian code
 */

#ifndef SCALARTYPEDEFS_HPP
#define SCALARTYPEDEFS_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "ScalarField/ScalarFieldTypedefs.hpp"

namespace GeoMHDiSCC {

   /// Namespace containing the required typedefs for a cartesian simulation
   namespace Datatypes {

      //
      // 1D case
      //

      /// Typedef for the 1D forward scalar field for the first transform
      typedef FlatDScalar1D FwdScalar1DType;

      /// Typedef for the 1D backward scalar field for the first transform
      typedef FlatDScalar1D BwdScalar1DType;

      //
      // 2D case
      //

      /// Typedef for the 2D forward scalar field for the first transform
      typedef FlatDScalar2D Fwd1DScalar2DType;

      /// Typedef for the 2D backward scalar field for the first transform
      typedef FlatDScalar2D Bwd1DScalar2DType;

      /// Typedef for the 2D forward scalar field for the second transform
      typedef FlatDScalar2D Fwd2DScalar2DType;

      /// Typedef for the 2D backward scalar field for the second transform
      typedef FlatDScalar2D Bwd2DScalar2DType;

      //
      // 3D case
      //

      /// Typedef for the 3D forward scalar field for the first transform
      typedef FlatZScalar3D Fwd1DScalar3DType;

      /// Typedef for the 3D backward scalar field for the first transform
      typedef FlatZScalar3D Bwd1DScalar3DType;

      /// Typedef for the 3D forward scalar field for the second transform
      typedef FlatDScalar3D Fwd2DScalar3DType;

      /// Typedef for the 3D backward scalar field for the second transform
      typedef FlatZScalar3D Bwd2DScalar3DType;

      /// Typedef for the 3D forward scalar field for the third transform
      typedef FlatDScalar3D Fwd3DScalar3DType;

      /// Typedef for the 3D backward scalar field for the third transform
      typedef FlatDScalar3D Bwd3DScalar3DType;

      // 1D case
      #if defined EPMPHOENIX_CODEDIM_1D

         /// Typedef for the physical space scalar
         typedef FwdScalar1DType PhysicalScalarType;

         /// Typedef for the spectral space scalar
         typedef BwdScalar1DType SpectralScalarType;

      // 2D case
      #elif defined EPMPHOENIX_CODEDIM_2D

         /// Typedef for the physical space scalar
         typedef Fwd2DScalar2DType PhysicalScalarType;

         /// Typedef for the spectral space scalar
         typedef Bwd1DScalar2DType SpectralScalarType;

      // 3D case
      #elif defined EPMPHOENIX_CODEDIM_3D

         /// Typedef for the physical space scalar
         typedef Fwd3DScalar3DType PhysicalScalarType;

         /// Typedef for the spectral space scalar
         typedef Bwd1DScalar3DType SpectralScalarType;

      // unspecified case
      #else
         #error "Number of dimension is simulation code not specified"
      #endif // EPMPHOENIX_CODEDIM_1D

   }
}

#endif // SCALARTYPEDEFS_HPP
