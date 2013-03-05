/** \file SpectralSelector.hpp
 *  \brief Small template to select the right spectral operator types
 *
 *  \mhdBug Needs test
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
#include "SpectralOperators/ChebyshevBoundary.hpp"
#include "SpectralOperators/UnitOperator.hpp"

namespace GeoMHDiSCC {

   namespace Spectral {

      template <Dimensions::Simulation::Id TId>  struct SpectralSelector;

      template <> struct SpectralSelector<Dimensions::Simulation::SIM1D>
      {
         /// Typedef for the spectral operator
         typedef  ChebyshevOperator  OpType;

         /// Typedef for the spectral boundary operator
         typedef  ChebyshevBoundary  BcType;
      };

      template <> struct SpectralSelector<Dimensions::Simulation::SIM2D>
      {
      };

      template <> struct SpectralSelector<Dimensions::Simulation::SIM3D>
      {
         /// Typedef for the spectral operator
         typedef  ChebyshevOperator  OpType;

         /// Typedef for the spectral boundary operator
         typedef  ChebyshevBoundary  BcType;
      };

   }
}

#endif // SPECTRALSELECTOR_HPP
