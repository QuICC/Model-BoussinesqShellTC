/** \file SpectralSelector.hpp
 *  \brief Small template to select the right spectral operator types
 */

#ifndef SPECTRALSELECTOR_HPP
#define SPECTRALSELECTOR_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "SpectralOperators/ChebyshevOperator.hpp"
#include "SpectralOperators/UnitOperator.hpp"

namespace GeoMHDiSCC {

   namespace Spectral {

      template <Dimensions::Transform::Id TId>  struct SpectralSelector;

      template <> struct SpectralSelector<Dimensions::Transform::TRA1D>
      {
         typedef  ChebyshevOperator  Type;
      };

      template <> struct SpectralSelector<Dimensions::Transform::TRA2D>
      {
         typedef  UnitOperator  Type;
      };

      template <> struct SpectralSelector<Dimensions::Transform::TRA3D>
      {
         typedef  ChebyshevOperator  Type;
      };

   }
}

#endif // SPECTRALSELECTOR_HPP
