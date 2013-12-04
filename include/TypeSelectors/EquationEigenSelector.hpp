/** 
 * @file EquationEigenSelector.hpp
 * @brief Selector to define equation eigen tools
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef EQUATIONEIGENSELECTOR_HPP
#define EQUATIONEIGENSELECTOR_HPP

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

         namespace EigenSelector = NoEigen;
      }
   }
#endif //GEOMHDISCC_SPATIALSCHEME_TTT

// Configure code to use TFT scheme
#ifdef GEOMHDISCC_SPATIALSCHEME_TFT
   
   #include "Equations/Tools/EquationEigen1DTools.hpp"

   namespace GeoMHDiSCC {

      namespace Equations {

         namespace EigenSelector = Eigen1D;
      }
   }
#endif //GEOMHDISCC_SPATIALSCHEME_TFT

// Configure code to use TFF scheme
#ifdef GEOMHDISCC_SPATIALSCHEME_TFF

   #include "Equations/Tools/EquationEigen2DTools.hpp"

   namespace GeoMHDiSCC {

      namespace Equations {

         namespace EigenSelector = Eigen2D;
      }
   }
#endif //GEOMHDISCC_SPATIALSCHEME_TFF

// Configure code to use FFF scheme
#ifdef GEOMHDISCC_SPATIALSCHEME_FFF

   #include "Equations/Tools/EquationEigen3DTools.hpp"

   namespace GeoMHDiSCC {

      namespace Equations {

         namespace EigenSelector = Eigen3D;
      }
   }
#endif //GEOMHDISCC_SPATIALSCHEME_FFF

// Configure code to use CFT scheme
#ifdef GEOMHDISCC_SPATIALSCHEME_CFT

   #include "Equations/Tools/EquationEigen1DTools.hpp"

   namespace GeoMHDiSCC {

      namespace Equations {

         namespace EigenSelector = Eigen1D;
      }
   }
#endif //GEOMHDISCC_SPATIALSCHEME_CFT

// Configure code to use AFT scheme
#ifdef GEOMHDISCC_SPATIALSCHEME_AFT

   #include "Equations/Tools/EquationEigen1DTools.hpp"

   namespace GeoMHDiSCC {

      namespace Equations {

         namespace EigenSelector = Eigen1D;
      }
   }
#endif //GEOMHDISCC_SPATIALSCHEME_AFT

// Configure code to use BLF scheme
#ifdef GEOMHDISCC_SPATIALSCHEME_BLF

   #include "Equations/Tools/EquationEigen2DTools.hpp"

   namespace GeoMHDiSCC {

      namespace Equations {

         namespace EigenSelector = Eigen2D;
      }
   }
#endif //GEOMHDISCC_SPATIALSCHEME_BLF

// Configure code to use SLF scheme
#ifdef GEOMHDISCC_SPATIALSCHEME_SLF

   #include "Equations/Tools/EquationEigen2DTools.hpp"

   namespace GeoMHDiSCC {

      namespace Equations {

         namespace EigenSelector = Eigen2D;
      }
   }
#endif //GEOMHDISCC_SPATIALSCHEME_SLF

// Configure code to use WFT scheme
#ifdef GEOMHDISCC_SPATIALSCHEME_WFT

   #include "Equations/Tools/EquationEigen1DTools.hpp"

   namespace GeoMHDiSCC {

      namespace Equations {

         namespace EigenSelector = Eigen1D;
      }
   }
#endif //GEOMHDISCC_SPATIALSCHEME_WFT

// Configure code to use WLF scheme
#ifdef GEOMHDISCC_SPATIALSCHEME_WLF

   #include "Equations/Tools/EquationEigen2DTools.hpp"

   namespace GeoMHDiSCC {

      namespace Equations {

         namespace EigenSelector = Eigen2D;
      }
   }
#endif //GEOMHDISCC_SPATIALSCHEME_WLF

#endif // EQUATIONEIGENSELECTOR_HPP
