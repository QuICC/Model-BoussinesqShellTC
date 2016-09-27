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
         typedef TTTScheme SpatialSelector;
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
         typedef TFTScheme SpatialSelector;
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
         typedef TFFScheme SpatialSelector;
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
         typedef FFFScheme SpatialSelector;
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
         typedef CFTScheme SpatialSelector;
      }
   }
#endif //GEOMHDISCC_SPATIALSCHEME_CFT

// Configure code to use AFT scheme
#ifdef GEOMHDISCC_SPATIALSCHEME_AFT
   // CFT includes
   #include "SpatialSchemes/3D/AFTScheme.hpp"

   namespace GeoMHDiSCC {

      namespace Schemes {

         /// Typedef for the CFT spatial scheme
         typedef AFTScheme SpatialSelector;
      }
   }
#endif //GEOMHDISCC_SPATIALSCHEME_AFT

// Configure code to use BLFl scheme
#ifdef GEOMHDISCC_SPATIALSCHEME_BLFL
   // BLFl includes
   #include "SpatialSchemes/3D/BLFlScheme.hpp"

   namespace GeoMHDiSCC {

      namespace Schemes {

         /// Typedef for the BLFl spatial scheme
         typedef BLFlScheme SpatialSelector;
      }
   }
#endif //GEOMHDISCC_SPATIALSCHEME_BLFL

// Configure code to use BLFm scheme
#ifdef GEOMHDISCC_SPATIALSCHEME_BLFM
   // BLFm includes
   #include "SpatialSchemes/3D/BLFmScheme.hpp"

   namespace GeoMHDiSCC {

      namespace Schemes {

         /// Typedef for the BLFm spatial scheme
         typedef BLFmScheme SpatialSelector;
      }
   }
#endif //GEOMHDISCC_SPATIALSCHEME_BLFM

// Configure code to use SLFl scheme
#ifdef GEOMHDISCC_SPATIALSCHEME_SLFL
   // SLFl includes
   #include "SpatialSchemes/3D/SLFlScheme.hpp"

   namespace GeoMHDiSCC {

      namespace Schemes {

         /// Typedef for the SLFl spatial scheme
         typedef SLFlScheme SpatialSelector;
      }
   }
#endif //GEOMHDISCC_SPATIALSCHEME_SLFL

// Configure code to use SLFm scheme
#ifdef GEOMHDISCC_SPATIALSCHEME_SLFM
   // SLFm includes
   #include "SpatialSchemes/3D/SLFmScheme.hpp"

   namespace GeoMHDiSCC {

      namespace Schemes {

         /// Typedef for the SLFm spatial scheme
         typedef SLFmScheme SpatialSelector;
      }
   }
#endif //GEOMHDISCC_SPATIALSCHEME_SLFM

// Configure code to use WFT scheme
#ifdef GEOMHDISCC_SPATIALSCHEME_WFT
   // WFT includes
   #include "SpatialSchemes/3D/WFTScheme.hpp"

   namespace GeoMHDiSCC {

      namespace Schemes {

         /// Typedef for the WFT spatial scheme
         typedef WFTScheme SpatialSelector;
      }
   }
#endif //GEOMHDISCC_SPATIALSCHEME_WFT

// Configure code to use WLFl scheme
#ifdef GEOMHDISCC_SPATIALSCHEME_WLFL
   // WLFl includes
   #include "SpatialSchemes/3D/WLFlScheme.hpp"

   namespace GeoMHDiSCC {

      namespace Schemes {

         /// Typedef for the WLFl spatial scheme
         typedef WLFlScheme SpatialSelector;
      }
   }
#endif //GEOMHDISCC_SPATIALSCHEME_WLFL

// Configure code to use WLFm scheme
#ifdef GEOMHDISCC_SPATIALSCHEME_WLFM
   // WLFm includes
   #include "SpatialSchemes/3D/WLFmScheme.hpp"

   namespace GeoMHDiSCC {

      namespace Schemes {

         /// Typedef for the WLFl spatial scheme
         typedef WLFmScheme SpatialSelector;
      }
   }
#endif //GEOMHDISCC_SPATIALSCHEME_WLFM

// Configure code to use TT scheme
#ifdef GEOMHDISCC_SPATIALSCHEME_TT
   
   // TTT includes
   #include "SpatialSchemes/2D/TTScheme.hpp"

   namespace GeoMHDiSCC {

      namespace Schemes {

         /// Typedef for the TT spatial scheme
         typedef TTScheme SpatialSelector;
      }
   }
#endif //GEOMHDISCC_SPATIALSCHEME_TT

// Configure code to use TF scheme
#ifdef GEOMHDISCC_SPATIALSCHEME_TF
   
   // TF includes
   #include "SpatialSchemes/2D/TFScheme.hpp"

   namespace GeoMHDiSCC {

      namespace Schemes {

         /// Typedef for the TF spatial scheme
         typedef TFScheme SpatialSelector;
      }
   }
#endif //GEOMHDISCC_SPATIALSCHEME_TF

#endif // SPATIALSCHEMESELECTOR_HPP
