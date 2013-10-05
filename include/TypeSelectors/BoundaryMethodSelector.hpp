/** 
 * @file BoundaryMethodSelector.hpp
 * @brief Small template to select the right spectral operator types
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef BOUNDARYMETHODSELECTOR_HPP
#define BOUNDARYMETHODSELECTOR_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//

namespace GeoMHDiSCC {

   namespace Boundary {

      template <Dimensions::Simulation::Id TId>  struct MethodSelector;
   }
}

   // Configure code to use TTT scheme
   #ifdef GEOMHDISCC_SPATIALSCHEME_TTT

      #include "SpectralOperators/TauChebyshev.hpp"
      #include "SpectralOperators/GalerkinChebyshev.hpp"

      namespace GeoMHDiSCC {
         namespace Boundary {

            template <> struct MethodSelector<Dimensions::Simulation::SIM1D>
            {
               #ifdef GEOMHDISCC_BOUNDARYMETHOD_GALERKIN
                  typedef  Spectral::GalerkinChebyshev  Type;
               #else
                  typedef  Spectral::TauChebyshev  Type;
               #endif //GEOMHDISCC_BOUNDARYMETHOD_GALERKIN
            };

            template <> struct MethodSelector<Dimensions::Simulation::SIM2D>
            {
               #ifdef GEOMHDISCC_BOUNDARYMETHOD_GALERKIN
                  typedef  Spectral::GalerkinChebyshev  Type;
               #else
                  typedef  Spectral::TauChebyshev  Type;
               #endif //GEOMHDISCC_BOUNDARYMETHOD_GALERKIN
            };

            template <> struct MethodSelector<Dimensions::Simulation::SIM3D>
            {
               /// Typedef for the spectral boundary operator
               #ifdef GEOMHDISCC_BOUNDARYMETHOD_GALERKIN
                  typedef  Spectral::GalerkinChebyshev  Type;
               #else
                  typedef  Spectral::TauChebyshev  Type;
               #endif //GEOMHDISCC_BOUNDARYMETHOD_GALERKIN
            };
         }
      }
   #endif //GEOMHDISCC_SPATIALSCHEME_TTT

   // Configure code to use TFT scheme
   #ifdef GEOMHDISCC_SPATIALSCHEME_TFT

      #include "SpectralOperators/TauChebyshev.hpp"
      #include "SpectralOperators/GalerkinChebyshev.hpp"

      namespace GeoMHDiSCC {
         namespace Boundary {

            template <> struct MethodSelector<Dimensions::Simulation::SIM1D>
            {
               #ifdef GEOMHDISCC_BOUNDARYMETHOD_GALERKIN
                  typedef  Spectral::GalerkinChebyshev  Type;
               #else
                  typedef  Spectral::TauChebyshev  Type;
               #endif //GEOMHDISCC_BOUNDARYMETHOD_GALERKIN
            };

            template <> struct MethodSelector<Dimensions::Simulation::SIM2D>
            {
            };

            template <> struct MethodSelector<Dimensions::Simulation::SIM3D>
            {
               #ifdef GEOMHDISCC_BOUNDARYMETHOD_GALERKIN
                  typedef  Spectral::GalerkinChebyshev  Type;
               #else
                  typedef  Spectral::TauChebyshev  Type;
               #endif //GEOMHDISCC_BOUNDARYMETHOD_GALERKIN
            };
         }
      }
   #endif //GEOMHDISCC_SPATIALSCHEME_TFT

   // Configure code to use TFF scheme
   #ifdef GEOMHDISCC_SPATIALSCHEME_TFF

      #include "SpectralOperators/TauChebyshev.hpp"
      #include "SpectralOperators/GalerkinChebyshev.hpp"

      namespace GeoMHDiSCC {
         namespace Boundary {

            template <> struct MethodSelector<Dimensions::Simulation::SIM1D>
            {
               #ifdef GEOMHDISCC_BOUNDARYMETHOD_GALERKIN
                  typedef  Spectral::GalerkinChebyshev  Type;
               #else
                  typedef  Spectral::TauChebyshev  Type;
               #endif //GEOMHDISCC_BOUNDARYMETHOD_GALERKIN
            };

            template <> struct MethodSelector<Dimensions::Simulation::SIM2D>
            {
            };

            template <> struct MethodSelector<Dimensions::Simulation::SIM3D>
            {
            };
         }
      }
   #endif //GEOMHDISCC_SPATIALSCHEME_TFF

   // Configure code to use FFF scheme
   #ifdef GEOMHDISCC_SPATIALSCHEME_FFF

      #include "SpectralOperators/TauChebyshev.hpp"
      #include "SpectralOperators/GalerkinChebyshev.hpp"

      namespace GeoMHDiSCC {
         namespace Boundary {

            template <> struct MethodSelector<Dimensions::Simulation::SIM1D>
            {
            };

            template <> struct MethodSelector<Dimensions::Simulation::SIM2D>
            {
            };

            template <> struct MethodSelector<Dimensions::Simulation::SIM3D>
            {
            };
         }
      }
   #endif //GEOMHDISCC_SPATIALSCHEME_FFF

   // Configure code to use CFT scheme
   #ifdef GEOMHDISCC_SPATIALSCHEME_CFT

      #include "SpectralOperators/TauCShellChebyshev.hpp"
      #include "SpectralOperators/GalerkinCShellChebyshev.hpp"

      namespace GeoMHDiSCC {
         namespace Boundary {

            template <> struct MethodSelector<Dimensions::Simulation::SIM1D>
            {
               #ifdef GEOMHDISCC_BOUNDARYMETHOD_GALERKIN
                  typedef  Spectral::GalerkinCShellChebyshev  Type;
               #else
                  typedef  Spectral::TauCShellChebyshev  Type;
               #endif //GEOMHDISCC_BOUNDARYMETHOD_GALERKIN
            };

            template <> struct MethodSelector<Dimensions::Simulation::SIM2D>
            {
            };

            template <> struct MethodSelector<Dimensions::Simulation::SIM3D>
            {
               #ifdef GEOMHDISCC_BOUNDARYMETHOD_GALERKIN
                  typedef  Spectral::GalerkinChebyshev  Type;
               #else
                  typedef  Spectral::TauChebyshev  Type;
               #endif //GEOMHDISCC_BOUNDARYMETHOD_GALERKIN
            };
         }
      }
   #endif //GEOMHDISCC_SPATIALSCHEME_CFT

   // Configure code to use SLF scheme
   #ifdef GEOMHDISCC_SPATIALSCHEME_SLF

      #include "SpectralOperators/TauSShellChebyshev.hpp"
      #include "SpectralOperators/GalerkinSShellChebyshev.hpp"

      namespace GeoMHDiSCC {
         namespace Boundary {

            template <> struct MethodSelector<Dimensions::Simulation::SIM1D>
            {
               #ifdef GEOMHDISCC_BOUNDARYMETHOD_GALERKIN
                  typedef  Spectral::GalerkinSShellChebyshev  Type;
               #else
                  typedef  Spectral::TauSShellChebyshev  Type;
               #endif //GEOMHDISCC_BOUNDARYMETHOD_GALERKIN
            };

            template <> struct MethodSelector<Dimensions::Simulation::SIM2D>
            {
            };

            template <> struct MethodSelector<Dimensions::Simulation::SIM3D>
            {
            };
         }
      }
   #endif //GEOMHDISCC_SPATIALSCHEME_SLF

   // Configure code to use WFT scheme
   #ifdef GEOMHDISCC_SPATIALSCHEME_WFT

      #include "SpectralOperators/TauCylindricalChebyshev.hpp"
      #include "SpectralOperators/GalerkinCylindricalChebyshev.hpp"

      namespace GeoMHDiSCC {
         namespace Boundary {

            template <> struct MethodSelector<Dimensions::Simulation::SIM1D>
            {
               #ifdef GEOMHDISCC_BOUNDARYMETHOD_GALERKIN
                  typedef  Spectral::GalerkinCylindricalChebyshev  Type;
               #else
                  typedef  Spectral::TauCylindricalChebyshev  Type;
               #endif //GEOMHDISCC_BOUNDARYMETHOD_GALERKIN
            };

            template <> struct MethodSelector<Dimensions::Simulation::SIM2D>
            {
            };

            template <> struct MethodSelector<Dimensions::Simulation::SIM3D>
            {
               #ifdef GEOMHDISCC_BOUNDARYMETHOD_GALERKIN
                  typedef  Spectral::GalerkinChebyshev  Type;
               #else
                  typedef  Spectral::TauChebyshev  Type;
               #endif //GEOMHDISCC_BOUNDARYMETHOD_GALERKIN
            };
         }
      }
   #endif //GEOMHDISCC_SPATIALSCHEME_WFT

   // Configure code to use WLF scheme
   #ifdef GEOMHDISCC_SPATIALSCHEME_WLF

      #include "SpectralOperators/TauSphericalChebyshev.hpp"
      #include "SpectralOperators/GalerkinSphericalChebyshev.hpp"

      namespace GeoMHDiSCC {
         namespace Boundary {

            template <> struct MethodSelector<Dimensions::Simulation::SIM1D>
            {
               #ifdef GEOMHDISCC_BOUNDARYMETHOD_GALERKIN
                  typedef  Spectral::GalerkinSphericalChebyshev  Type;
               #else
                  typedef  Spectral::TauSphericalChebyshev  Type;
               #endif //GEOMHDISCC_BOUNDARYMETHOD_GALERKIN
            };

            template <> struct MethodSelector<Dimensions::Simulation::SIM2D>
            {
            };

            template <> struct MethodSelector<Dimensions::Simulation::SIM3D>
            {
            };
         }
      }
   #endif //GEOMHDISCC_SPATIALSCHEME_WLF

#endif // BOUNDARYMETHODSELECTOR_HPP
