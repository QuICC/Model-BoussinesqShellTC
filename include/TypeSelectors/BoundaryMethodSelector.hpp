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
#include <tr1/tuple>

// External includes
//

// Project includes
//
#include "BoundaryCondition/BoundaryCoordinator.hpp"

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

         typedef MHDFloat  BCIndex;

         typedef BoundaryCoordinator<BCIndex, MethodSelector<Dimensions::Simulation::SIM1D>::Type,MethodSelector<Dimensions::Simulation::SIM2D>::Type,MethodSelector<Dimensions::Simulation::SIM3D>::Type> CoordinatorSelector;

         /// Flag to specify index independent boundary conditions
         const BCIndex INDEPENDENT = std::numeric_limits<MHDFloat>::min();
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

         typedef MHDFloat   BCIndex;

         typedef BoundaryCoordinator<BCIndex, MethodSelector<Dimensions::Simulation::SIM1D>::Type, int, MethodSelector<Dimensions::Simulation::SIM3D>::Type> CoordinatorSelector;

         /// Flag to specify index independent boundary conditions
         const BCIndex INDEPENDENT = std::numeric_limits<MHDFloat>::min();
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

         typedef std::tr1::tuple<MHDFloat,MHDFloat>   BCIndex;

         typedef BoundaryCoordinator<BCIndex, MethodSelector<Dimensions::Simulation::SIM1D>::Type, int, int> CoordinatorSelector;

         /// Flag to specify index independent boundary conditions
         const BCIndex INDEPENDENT;
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

         typedef std::tr1::tuple<MHDFloat,MHDFloat,MHDFloat>   BCIndex;

         typedef BoundaryCoordinator<BCIndex, int, int, int> CoordinatorSelector;

         /// Flag to specify index independent boundary conditions
         const BCIndex INDEPENDENT;
      }
   }
#endif //GEOMHDISCC_SPATIALSCHEME_FFF

// Configure code to use CFT scheme
#ifdef GEOMHDISCC_SPATIALSCHEME_CFT

   #include "SpectralOperators/TauChebyshev.hpp"
   #include "SpectralOperators/GalerkinChebyshev.hpp"
   #include "SpectralOperators/TauCylinderChebyshev.hpp"
   #include "SpectralOperators/GalerkinCylinderChebyshev.hpp"

   namespace GeoMHDiSCC {

      namespace Boundary {

         template <> struct MethodSelector<Dimensions::Simulation::SIM1D>
         {
            #ifdef GEOMHDISCC_BOUNDARYMETHOD_GALERKIN
               typedef  Spectral::GalerkinCylinderChebyshev  Type;
            #else
               typedef  Spectral::TauCylinderChebyshev  Type;
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

         typedef MHDFloat   BCIndex;

         typedef BoundaryCoordinator<BCIndex, MethodSelector<Dimensions::Simulation::SIM1D>::Type, int, MethodSelector<Dimensions::Simulation::SIM3D>::Type> CoordinatorSelector;

         /// Flag to specify index independent boundary conditions
         const BCIndex INDEPENDENT = std::numeric_limits<MHDFloat>::min();
      }
   }
#endif //GEOMHDISCC_SPATIALSCHEME_CFT

// Configure code to use AFT scheme
#ifdef GEOMHDISCC_SPATIALSCHEME_AFT

   #include "SpectralOperators/TauChebyshev.hpp"
   #include "SpectralOperators/GalerkinChebyshev.hpp"
   #include "SpectralOperators/TauAnnulusChebyshev.hpp"
   #include "SpectralOperators/GalerkinAnnulusChebyshev.hpp"

   namespace GeoMHDiSCC {

      namespace Boundary {

         template <> struct MethodSelector<Dimensions::Simulation::SIM1D>
         {
            #ifdef GEOMHDISCC_BOUNDARYMETHOD_GALERKIN
               typedef  Spectral::GalerkinAnnulusChebyshev  Type;
            #else
               typedef  Spectral::TauAnnulusChebyshev  Type;
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

         typedef MHDFloat   BCIndex;

         typedef BoundaryCoordinator<BCIndex, MethodSelector<Dimensions::Simulation::SIM1D>::Type, int, MethodSelector<Dimensions::Simulation::SIM3D>::Type> CoordinatorSelector;

         /// Flag to specify index independent boundary conditions
         const BCIndex INDEPENDENT = std::numeric_limits<MHDFloat>::min();
      }
   }
#endif //GEOMHDISCC_SPATIALSCHEME_AFT

// Configure code to use BLF scheme
#ifdef GEOMHDISCC_SPATIALSCHEME_BLF

   #include "SpectralOperators/TauSphereChebyshev.hpp"
   #include "SpectralOperators/GalerkinSphereChebyshev.hpp"

   namespace GeoMHDiSCC {

      namespace Boundary {

         template <> struct MethodSelector<Dimensions::Simulation::SIM1D>
         {
            #ifdef GEOMHDISCC_BOUNDARYMETHOD_GALERKIN
               typedef  Spectral::GalerkinSphereChebyshev  Type;
            #else
               typedef  Spectral::TauSphereChebyshev  Type;
            #endif //GEOMHDISCC_BOUNDARYMETHOD_GALERKIN
         };

         template <> struct MethodSelector<Dimensions::Simulation::SIM2D>
         {
         };

         template <> struct MethodSelector<Dimensions::Simulation::SIM3D>
         {
         };

         typedef std::tr1::tuple<MHDFloat,MHDFloat>   BCIndex;

         typedef BoundaryCoordinator<BCIndex, MethodSelector<Dimensions::Simulation::SIM1D>::Type, int, int> CoordinatorSelector;

         /// Flag to specify index independent boundary conditions
         const BCIndex INDEPENDENT;
      }
   }
#endif //GEOMHDISCC_SPATIALSCHEME_BLF

// Configure code to use SLF scheme
#ifdef GEOMHDISCC_SPATIALSCHEME_SLF

   #include "SpectralOperators/TauShellChebyshev.hpp"
   #include "SpectralOperators/GalerkinShellChebyshev.hpp"

   namespace GeoMHDiSCC {

      namespace Boundary {

         template <> struct MethodSelector<Dimensions::Simulation::SIM1D>
         {
            #ifdef GEOMHDISCC_BOUNDARYMETHOD_GALERKIN
               typedef  Spectral::GalerkinShellChebyshev  Type;
            #else
               typedef  Spectral::TauShellChebyshev  Type;
            #endif //GEOMHDISCC_BOUNDARYMETHOD_GALERKIN
         };

         template <> struct MethodSelector<Dimensions::Simulation::SIM2D>
         {
         };

         template <> struct MethodSelector<Dimensions::Simulation::SIM3D>
         {
         };

         typedef std::tr1::tuple<MHDFloat,MHDFloat>   BCIndex;

         typedef BoundaryCoordinator<BCIndex, MethodSelector<Dimensions::Simulation::SIM1D>::Type, int, int> CoordinatorSelector;

         /// Flag to specify index independent boundary conditions
         const BCIndex INDEPENDENT;
      }
   }
#endif //GEOMHDISCC_SPATIALSCHEME_SLF

// Configure code to use WFT scheme
#ifdef GEOMHDISCC_SPATIALSCHEME_WFT

   #include "SpectralOperators/TauCylinderWorland.hpp"
   #include "SpectralOperators/GalerkinCylinderWorland.hpp"

   namespace GeoMHDiSCC {

      namespace Boundary {

         template <> struct MethodSelector<Dimensions::Simulation::SIM1D>
         {
            #ifdef GEOMHDISCC_BOUNDARYMETHOD_GALERKIN
               typedef  Spectral::GalerkinCylinderWorland  Type;
            #else
               typedef  Spectral::TauCylinderWorland  Type;
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

   #include "SpectralOperators/TauSphereWorland.hpp"
   #include "SpectralOperators/GalerkinSphereWorland.hpp"

   namespace GeoMHDiSCC {

      namespace Boundary {

         template <> struct MethodSelector<Dimensions::Simulation::SIM1D>
         {
            #ifdef GEOMHDISCC_BOUNDARYMETHOD_GALERKIN
               typedef  Spectral::GalerkinSphereWorland  Type;
            #else
               typedef  Spectral::TauSphereWorland  Type;
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
