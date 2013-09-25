/** 
 * @file EquationToolsSelector.hpp
 * @brief Selector to define equation tools typedef
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef EQUATIONTOOLSSELECTOR_HPP
#define EQUATIONTOOLSSELECTOR_HPP

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
   
   #include "Equations/Tools/EquationNoEigenTools.hpp"

   namespace GeoMHDiSCC {

      namespace Equations {

         /// Typedef for the TTT spatial scheme
         typedef EquationNoEigenTools EquationToolsType;
      }
   }
#endif //GEOMHDISCC_SPATIALSCHEME_TTT

// Configure code to use TFT scheme
#ifdef GEOMHDISCC_SPATIALSCHEME_TFT
   
   #include "Equations/Tools/Equation1DEigenTools.hpp"

   namespace GeoMHDiSCC {

      namespace Equations {

         /// Typedef for the TTT spatial scheme
         typedef Equation1DEigenTools EquationToolsType;
      }
   }
#endif //GEOMHDISCC_SPATIALSCHEME_TFT

// Configure code to use TFF scheme
#ifdef GEOMHDISCC_SPATIALSCHEME_TFF

   #include "Equations/Tools/Equation2DEigenTools.hpp"

   namespace GeoMHDiSCC {

      namespace Equations {

         /// Typedef for the TTT spatial scheme
         typedef Equation2DEigenTools EquationToolsType;
      }
   }
#endif //GEOMHDISCC_SPATIALSCHEME_TFF

// Configure code to use FFF scheme
#ifdef GEOMHDISCC_SPATIALSCHEME_FFF

   #include "Equations/Tools/Equation3DEigenTools.hpp"

   namespace GeoMHDiSCC {

      namespace Equations {

         /// Typedef for the TTT spatial scheme
         typedef Equation3DEigenTools EquationToolsType;
      }
   }
#endif //GEOMHDISCC_SPATIALSCHEME_FFF

// Configure code to use CFT scheme
#ifdef GEOMHDISCC_SPATIALSCHEME_CFT

   #include "Equations/Tools/Equation1DEigenTools.hpp"

   namespace GeoMHDiSCC {

      namespace Equations {

         /// Typedef for the TTT spatial scheme
         typedef Equation1DEigenTools EquationToolsType;
      }
   }
#endif //GEOMHDISCC_SPATIALSCHEME_CFT

// Configure code to use SLF scheme
#ifdef GEOMHDISCC_SPATIALSCHEME_SLF

   #include "Equations/Tools/Equation2DEigenTools.hpp"

   namespace GeoMHDiSCC {

      namespace Equations {

         /// Typedef for the TTT spatial scheme
         typedef Equation2DEigenTools EquationToolsType;
      }
   }
#endif //GEOMHDISCC_SPATIALSCHEME_SLF

// Configure code to use WFT scheme
#ifdef GEOMHDISCC_SPATIALSCHEME_WFT

   #include "Equations/Tools/Equation1DEigenTools.hpp"

   namespace GeoMHDiSCC {

      namespace Equations {

         /// Typedef for the TTT spatial scheme
         typedef Equation2DEigenTools EquationToolsType;
      }
   }
#endif //GEOMHDISCC_SPATIALSCHEME_WFT

// Configure code to use WLF scheme
#ifdef GEOMHDISCC_SPATIALSCHEME_WLF

   #include "Equations/Tools/Equation2DEigenTools.hpp"

   namespace GeoMHDiSCC {

      namespace Equations {

         /// Typedef for the TTT spatial scheme
         typedef Equation2DEigenTools EquationToolsType;
      }
   }
#endif //GEOMHDISCC_SPATIALSCHEME_WLF

#endif // EQUATIONTOOLSSELECTOR_HPP
