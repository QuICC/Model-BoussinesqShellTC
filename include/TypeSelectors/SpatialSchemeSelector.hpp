/** 
 * @file SpatialSchemeSelector.hpp
 * @brief Selector to define spatial scheme typedef
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef SPATIALSCHEMESELECTOR_HPP
#define SPATIALSCHEMESELECTOR_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//

// Configure code to use TTT scheme
#ifdef GEOMHDISCC_SPATIALSCHEME_TTT
   
   // TTT includes
   #include "SpatialSchemes/3D/TTTScheme.hpp"

   namespace GeoMHDiSCC {

      namespace Schemes {

         /// Typedef for the TTT spatial scheme
         typedef TTTScheme SpatialType;
      }
   }
#endif //GEOMHDISCC_SPATIALSCHEME_TTT

// Configure code to use TFT scheme
#ifdef GEOMHDISCC_SPATIALSCHEME_TFT
   
   // TFT includes
   #include "SpatialSchemes/3D/TFTScheme.hpp"

   namespace GeoMHDiSCC {

      namespace Schemes {

         /// Typedef for the TFT spatial scheme
         typedef TFTScheme SpatialType;
      }
   }
#endif //GEOMHDISCC_SPATIALSCHEME_TFT

// Configure code to use TFF scheme
#ifdef GEOMHDISCC_SPATIALSCHEME_TFF
   // TFF includes
   #include "SpatialSchemes/3D/TFFScheme.hpp"

   namespace GeoMHDiSCC {

      namespace Schemes {

         /// Typedef for the TFF spatial scheme
         typedef TFFScheme SpatialType;
      }
   }
#endif //GEOMHDISCC_SPATIALSCHEME_TFF

// Configure code to use FFF scheme
#ifdef GEOMHDISCC_SPATIALSCHEME_FFF
   // FFF includes
   #include "SpatialSchemes/3D/FFFScheme.hpp"

   namespace GeoMHDiSCC {

      namespace Schemes {

         /// Typedef for the FFF spatial scheme
         typedef FFFScheme SpatialType;
      }
   }
#endif //GEOMHDISCC_SPATIALSCHEME_FFF

// Configure code to use CFT scheme
#ifdef GEOMHDISCC_SPATIALSCHEME_CFT
   // CFT includes
   #include "SpatialSchemes/3D/CFTScheme.hpp"

   namespace GeoMHDiSCC {

      namespace Schemes {

         /// Typedef for the CFT spatial scheme
         typedef CFTScheme SpatialType;
      }
   }
#endif //GEOMHDISCC_SPATIALSCHEME_CFT

// Configure code to use WFT scheme
#ifdef GEOMHDISCC_SPATIALSCHEME_WFT
   // WFT includes
   #include "SpatialSchemes/3D/WFTScheme.hpp"

   namespace GeoMHDiSCC {

      namespace Schemes {

         /// Typedef for the WFT spatial scheme
         typedef WFTScheme SpatialType;
      }
   }
#endif //GEOMHDISCC_SPATIALSCHEME_WFT

// Configure code to use SLF scheme
#ifdef GEOMHDISCC_SPATIALSCHEME_SLF
   // SLF includes
   #include "SpatialSchemes/3D/SLFScheme.hpp"

   namespace GeoMHDiSCC {

      namespace Schemes {

         /// Typedef for the SLF spatial scheme
         typedef SLFScheme SpatialType;
      }
   }
#endif //GEOMHDISCC_SPATIALSCHEME_SLF

// Configure code to use WLF scheme
#ifdef GEOMHDISCC_SPATIALSCHEME_WLF
   // WLF includes
   #include "SpatialSchemes/3D/WLFScheme.hpp"

   namespace GeoMHDiSCC {

      namespace Schemes {

         /// Typedef for the WLF spatial scheme
         typedef WLFScheme SpatialType;
      }
   }
#endif //GEOMHDISCC_SPATIALSCHEME_WLF

#endif // SPATIALSCHEMESELECTOR_HPP
