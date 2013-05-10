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
#include "SpectralOperators/ChebyshevBoundary.hpp"
#include "SpectralOperators/UnitOperator.hpp"

namespace GeoMHDiSCC {

   namespace Spectral {

      template <Dimensions::Simulation::Id TId>  struct SpectralSelector;

      // Configure code to use TFT scheme
      #ifdef GEOMHDISCC_SPATIALSCHEME_TFT
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
      #endif //GEOMHDISCC_SPATIALSCHEME_TFT

      // Configure code to use TFF scheme
      #ifdef GEOMHDISCC_SPATIALSCHEME_TFF
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
         };
      #endif //GEOMHDISCC_SPATIALSCHEME_TFF

   }
}

#endif // SPECTRALSELECTOR_HPP
