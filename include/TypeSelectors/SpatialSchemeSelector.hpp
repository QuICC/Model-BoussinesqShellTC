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
#ifdef QUICC_SPATIALSCHEME_TTT
   
   // TTT includes
   #include "SpatialSchemes/3D/TTTScheme.hpp"

   namespace QuICC {

      namespace Schemes {

         /// Typedef for the TTT spatial scheme
         typedef TTTScheme SpatialSelector;
      }
   }
#endif //QUICC_SPATIALSCHEME_TTT

// Configure code to use TFT scheme
#ifdef QUICC_SPATIALSCHEME_TFT
   
   // TFT includes
   #include "SpatialSchemes/3D/TFTScheme.hpp"

   namespace QuICC {

      namespace Schemes {

         /// Typedef for the TFT spatial scheme
         typedef TFTScheme SpatialSelector;
      }
   }
#endif //QUICC_SPATIALSCHEME_TFT

// Configure code to use TFF scheme
#ifdef QUICC_SPATIALSCHEME_TFF
   // TFF includes
   #include "SpatialSchemes/3D/TFFScheme.hpp"

   namespace QuICC {

      namespace Schemes {

         /// Typedef for the TFF spatial scheme
         typedef TFFScheme SpatialSelector;
      }
   }
#endif //QUICC_SPATIALSCHEME_TFF

// Configure code to use FFF scheme
#ifdef QUICC_SPATIALSCHEME_FFF
   // FFF includes
   #include "SpatialSchemes/3D/FFFScheme.hpp"

   namespace QuICC {

      namespace Schemes {

         /// Typedef for the FFF spatial scheme
         typedef FFFScheme SpatialSelector;
      }
   }
#endif //QUICC_SPATIALSCHEME_FFF

// Configure code to use CFT scheme
#ifdef QUICC_SPATIALSCHEME_CFT
   // CFT includes
   #include "SpatialSchemes/3D/CFTScheme.hpp"

   namespace QuICC {

      namespace Schemes {

         /// Typedef for the CFT spatial scheme
         typedef CFTScheme SpatialSelector;
      }
   }
#endif //QUICC_SPATIALSCHEME_CFT

// Configure code to use AFT scheme
#ifdef QUICC_SPATIALSCHEME_AFT
   // CFT includes
   #include "SpatialSchemes/3D/AFTScheme.hpp"

   namespace QuICC {

      namespace Schemes {

         /// Typedef for the CFT spatial scheme
         typedef AFTScheme SpatialSelector;
      }
   }
#endif //QUICC_SPATIALSCHEME_AFT

// Configure code to use BLFl scheme
#ifdef QUICC_SPATIALSCHEME_BLFL
   // BLFl includes
   #include "SpatialSchemes/3D/BLFlScheme.hpp"

   namespace QuICC {

      namespace Schemes {

         /// Typedef for the BLFl spatial scheme
         typedef BLFlScheme SpatialSelector;
      }
   }
#endif //QUICC_SPATIALSCHEME_BLFL

// Configure code to use BLFm scheme
#ifdef QUICC_SPATIALSCHEME_BLFM
   // BLFm includes
   #include "SpatialSchemes/3D/BLFmScheme.hpp"

   namespace QuICC {

      namespace Schemes {

         /// Typedef for the BLFm spatial scheme
         typedef BLFmScheme SpatialSelector;
      }
   }
#endif //QUICC_SPATIALSCHEME_BLFM

// Configure code to use SLFl scheme
#ifdef QUICC_SPATIALSCHEME_SLFL
   // SLFl includes
   #include "SpatialSchemes/3D/SLFlScheme.hpp"

   namespace QuICC {

      namespace Schemes {

         /// Typedef for the SLFl spatial scheme
         typedef SLFlScheme SpatialSelector;
      }
   }
#endif //QUICC_SPATIALSCHEME_SLFL

// Configure code to use SLFm scheme
#ifdef QUICC_SPATIALSCHEME_SLFM
   // SLFm includes
   #include "SpatialSchemes/3D/SLFmScheme.hpp"

   namespace QuICC {

      namespace Schemes {

         /// Typedef for the SLFm spatial scheme
         typedef SLFmScheme SpatialSelector;
      }
   }
#endif //QUICC_SPATIALSCHEME_SLFM

// Configure code to use WFT scheme
#ifdef QUICC_SPATIALSCHEME_WFT
   // WFT includes
   #include "SpatialSchemes/3D/WFTScheme.hpp"

   namespace QuICC {

      namespace Schemes {

         /// Typedef for the WFT spatial scheme
         typedef WFTScheme SpatialSelector;
      }
   }
#endif //QUICC_SPATIALSCHEME_WFT

// Configure code to use WLFl scheme
#ifdef QUICC_SPATIALSCHEME_WLFL
   // WLFl includes
   #include "SpatialSchemes/3D/WLFlScheme.hpp"

   namespace QuICC {

      namespace Schemes {

         /// Typedef for the WLFl spatial scheme
         typedef WLFlScheme SpatialSelector;
      }
   }
#endif //QUICC_SPATIALSCHEME_WLFL

// Configure code to use WLFm scheme
#ifdef QUICC_SPATIALSCHEME_WLFM
   // WLFm includes
   #include "SpatialSchemes/3D/WLFmScheme.hpp"

   namespace QuICC {

      namespace Schemes {

         /// Typedef for the WLFl spatial scheme
         typedef WLFmScheme SpatialSelector;
      }
   }
#endif //QUICC_SPATIALSCHEME_WLFM

// Configure code to use TT scheme
#ifdef QUICC_SPATIALSCHEME_TT
   
   // TTT includes
   #include "SpatialSchemes/2D/TTScheme.hpp"

   namespace QuICC {

      namespace Schemes {

         /// Typedef for the TT spatial scheme
         typedef TTScheme SpatialSelector;
      }
   }
#endif //QUICC_SPATIALSCHEME_TT

// Configure code to use TF scheme
#ifdef QUICC_SPATIALSCHEME_TF
   
   // TF includes
   #include "SpatialSchemes/2D/TFScheme.hpp"

   namespace QuICC {

      namespace Schemes {

         /// Typedef for the TF spatial scheme
         typedef TFScheme SpatialSelector;
      }
   }
#endif //QUICC_SPATIALSCHEME_TF

#endif // SPATIALSCHEMESELECTOR_HPP
