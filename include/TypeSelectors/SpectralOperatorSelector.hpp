/** 
 * @file OperatorSelector.hpp
 * @brief Small template to select the right spectral operator types
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef SPECTRALOPERATORSELECTOR_HPP
#define SPECTRALOPERATORSELECTOR_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//

namespace GeoMHDiSCC {

   namespace Spectral {

      template <Dimensions::Simulation::Id TId>  struct OperatorSelector;
   }
}

   // Configure code to use TTT scheme
   #ifdef GEOMHDISCC_SPATIALSCHEME_TTT

      #include "SpectralOperators/ChebyshevOperator.hpp"

      namespace GeoMHDiSCC {
         namespace Spectral {

            template <> struct OperatorSelector<Dimensions::Simulation::SIM1D>
            {
               typedef  ChebyshevOperator  Type;
            };

            template <> struct OperatorSelector<Dimensions::Simulation::SIM2D>
            {
               typedef  ChebyshevOperator  Type;
            };

            template <> struct OperatorSelector<Dimensions::Simulation::SIM3D>
            {
               typedef  ChebyshevOperator  Type;
            };
         }
      }
   #endif //GEOMHDISCC_SPATIALSCHEME_TTT

   // Configure code to use TFT scheme
   #ifdef GEOMHDISCC_SPATIALSCHEME_TFT

      #include "SpectralOperators/ChebyshevOperator.hpp"

      namespace GeoMHDiSCC {
         namespace Spectral {

            template <> struct OperatorSelector<Dimensions::Simulation::SIM1D>
            {
               typedef  ChebyshevOperator  Type;
            };

            template <> struct OperatorSelector<Dimensions::Simulation::SIM2D>
            {
            };

            template <> struct OperatorSelector<Dimensions::Simulation::SIM3D>
            {
               typedef  ChebyshevOperator  Type;
            };
         }
      }
   #endif //GEOMHDISCC_SPATIALSCHEME_TFT

   // Configure code to use TFF scheme
   #ifdef GEOMHDISCC_SPATIALSCHEME_TFF

      #include "SpectralOperators/ChebyshevOperator.hpp"

      namespace GeoMHDiSCC {
         namespace Spectral {

            template <> struct OperatorSelector<Dimensions::Simulation::SIM1D>
            {
               /// Typedef for the spectral operator
               typedef  ChebyshevOperator  Type;
            };

            template <> struct OperatorSelector<Dimensions::Simulation::SIM2D>
            {
            };

            template <> struct OperatorSelector<Dimensions::Simulation::SIM3D>
            {
            };
         }
      }
   #endif //GEOMHDISCC_SPATIALSCHEME_TFF

   // Configure code to use FFF scheme
   #ifdef GEOMHDISCC_SPATIALSCHEME_FFF

      #include "SpectralOperators/ChebyshevOperator.hpp"

      namespace GeoMHDiSCC {
         namespace Spectral {

            template <> struct OperatorSelector<Dimensions::Simulation::SIM1D>
            {
            };

            template <> struct OperatorSelector<Dimensions::Simulation::SIM2D>
            {
            };

            template <> struct OperatorSelector<Dimensions::Simulation::SIM3D>
            {
            };
         }
      }
   #endif //GEOMHDISCC_SPATIALSCHEME_FFF

   // Configure code to use CFT scheme
   #ifdef GEOMHDISCC_SPATIALSCHEME_CFT

      #include "SpectralOperators/ChebyshevOperator.hpp"
      #include "SpectralOperators/CShellChebyshevOperator.hpp"

      namespace GeoMHDiSCC {
         namespace Spectral {

            template <> struct OperatorSelector<Dimensions::Simulation::SIM1D>
            {
               typedef  CShellChebyshevOperator  Type;
            };

            template <> struct OperatorSelector<Dimensions::Simulation::SIM2D>
            {
            };

            template <> struct OperatorSelector<Dimensions::Simulation::SIM3D>
            {
               typedef  ChebyshevOperator  Type;
            };
         }
      }
   #endif //GEOMHDISCC_SPATIALSCHEME_CFT

   // Configure code to use SLF scheme
   #ifdef GEOMHDISCC_SPATIALSCHEME_SLF

      #include "SpectralOperators/SShellChebyshevOperator.hpp"

      namespace GeoMHDiSCC {
         namespace Spectral {

            template <> struct OperatorSelector<Dimensions::Simulation::SIM1D>
            {
               typedef  SShellChebyshevOperator  Type;
            };

            template <> struct OperatorSelector<Dimensions::Simulation::SIM2D>
            {
            };

            template <> struct OperatorSelector<Dimensions::Simulation::SIM3D>
            {
            };
         }
      }
   #endif //GEOMHDISCC_SPATIALSCHEME_SLF

   // Configure code to use WFT scheme
   #ifdef GEOMHDISCC_SPATIALSCHEME_WFT

      #include "SpectralOperators/ChebyshevOperator.hpp"
      #include "SpectralOperators/CylindricalChebyshevOperator.hpp"

      namespace GeoMHDiSCC {
         namespace Spectral {

            template <> struct OperatorSelector<Dimensions::Simulation::SIM1D>
            {
               typedef  CylindricalChebyshevOperator  Type;
            };

            template <> struct OperatorSelector<Dimensions::Simulation::SIM2D>
            {
            };

            template <> struct OperatorSelector<Dimensions::Simulation::SIM3D>
            {
               typedef  ChebyshevOperator  Type;
            };
         }
      }
   #endif //GEOMHDISCC_SPATIALSCHEME_WFT

   // Configure code to use WLF scheme
   #ifdef GEOMHDISCC_SPATIALSCHEME_WLF

      #include "SpectralOperators/SphericalChebyshevOperator.hpp"

      namespace GeoMHDiSCC {
         namespace Spectral {

            template <> struct OperatorSelector<Dimensions::Simulation::SIM1D>
            {
               typedef  SphericalChebyshevOperator  Type;
            };

            template <> struct OperatorSelector<Dimensions::Simulation::SIM2D>
            {
            };

            template <> struct OperatorSelector<Dimensions::Simulation::SIM3D>
            {
            };
         }
      }
   #endif //GEOMHDISCC_SPATIALSCHEME_WLF

#endif // SPECTRALOPERATORSELECTOR_HPP
