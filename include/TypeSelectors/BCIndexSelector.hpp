/** 
 * @file BCIndexSelector.hpp
 * @brief Selector to define equation tools typedef
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef BCINDEXSELECTOR_HPP
#define BCINDEXSELECTOR_HPP

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
   
   namespace GeoMHDiSCC {

      namespace Boundary {
         /// Typedef for an boundary condition index type
         typedef int BCIndex;
      }
   }
#endif //GEOMHDISCC_SPATIALSCHEME_TTT

// Configure code to use TFT scheme
#ifdef GEOMHDISCC_SPATIALSCHEME_TFT
   
   namespace GeoMHDiSCC {

      namespace Boundary {
         /// Typedef for an boundary condition index type
         typedef int BCIndex;
      }
   }
#endif //GEOMHDISCC_SPATIALSCHEME_TFT

// Configure code to use TFF scheme
#ifdef GEOMHDISCC_SPATIALSCHEME_TFF

   namespace GeoMHDiSCC {

      namespace Boundary {
         /// Typedef for an boundary condition index type
         typedef std::pair<int,int> BCIndex;
      }
   }
#endif //GEOMHDISCC_SPATIALSCHEME_TFF

// Configure code to use FFF scheme
#ifdef GEOMHDISCC_SPATIALSCHEME_FFF

   namespace GeoMHDiSCC {

      namespace Boundary {
         /// Typedef for an boundary condition index type
         typedef std::tr1::tuple<int,int,int> BCIndex;
      }
   }
#endif //GEOMHDISCC_SPATIALSCHEME_FFF

// Configure code to use CFT scheme
#ifdef GEOMHDISCC_SPATIALSCHEME_CFT

   namespace GeoMHDiSCC {

      namespace Boundary {
         /// Typedef for an boundary condition index type
         typedef int BCIndex;
      }
   }
#endif //GEOMHDISCC_SPATIALSCHEME_CFT

// Configure code to use SLF scheme
#ifdef GEOMHDISCC_SPATIALSCHEME_SLF

   namespace GeoMHDiSCC {

      namespace Boundary {
         /// Typedef for an boundary condition index type
         typedef std::pair<int,int> BCIndex;
      }
   }
#endif //GEOMHDISCC_SPATIALSCHEME_SLF

// Configure code to use WFT scheme
#ifdef GEOMHDISCC_SPATIALSCHEME_WFT

   namespace GeoMHDiSCC {

      namespace Boundary {
         /// Typedef for an boundary condition index type
         typedef int BCIndex;
      }
   }
#endif //GEOMHDISCC_SPATIALSCHEME_WFT

// Configure code to use WLF scheme
#ifdef GEOMHDISCC_SPATIALSCHEME_WLF

   namespace GeoMHDiSCC {

      namespace Boundary {
         /// Typedef for an boundary condition index type
         typedef std::pair<int,int> BCIndex;
      }
   }
#endif //GEOMHDISCC_SPATIALSCHEME_WLF

#endif // BCINDEXSELECTOR_HPP
