/** 
 * @file EquationEigenSelector.cpp
 * @brief Definitions of the eigen direction tools selector
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// System includes
//

// External includes
//

// Class include
//
#include "TypeSelectors/EquationEigenSelector.hpp"

// Project includes
//


// Configure code to use TTT scheme
#ifdef QUICC_SPATIALSCHEME_TTT
   
   #include "Equations/Tools/NoEigenTools.hpp"

   namespace GeoMHDiSCC {

      namespace Equations {

         SharedIEigenTools eigenSelector(const CouplingInformation::IndexType indexType)
         {
            if(indexType == CouplingInformation::SINGLE)
            {
               SharedNoEigenTools spTools(new NoEigenTools());

               return spTools;

            } else
            {
               throw Exception("Unknown eigen tools type");
            }
         }
      }
   }
#endif //QUICC_SPATIALSCHEME_TTT

// Configure code to use TFT scheme
#ifdef QUICC_SPATIALSCHEME_TFT
   
   #include "Equations/Tools/Eigen1DTools.hpp"

   namespace GeoMHDiSCC {

      namespace Equations {

         SharedIEigenTools eigenSelector(const CouplingInformation::IndexType indexType)
         {
            if(indexType == CouplingInformation::SLOWEST_SINGLE_RHS)
            {
               SharedEigen1DTools spTools(new Eigen1DTools());

               return spTools;

            } else
            {
               throw Exception("Unknown eigen tools type");
            }
         }

      }
   }
#endif //QUICC_SPATIALSCHEME_TFT

// Configure code to use TFF scheme
#ifdef QUICC_SPATIALSCHEME_TFF

   #include "Equations/Tools/Eigen2DTools.hpp"

   namespace GeoMHDiSCC {

      namespace Equations {

         SharedIEigenTools eigenSelector(const CouplingInformation::IndexType indexType)
         {
            if(indexType == CouplingInformation::MODE)
            {
               SharedEigen2DTools spTools(new Eigen2DTools());

               return spTools;

            } else
            {
               throw Exception("Unknown eigen tools type");
            }
         }

      }
   }
#endif //QUICC_SPATIALSCHEME_TFF

// Configure code to use FFF scheme
#ifdef QUICC_SPATIALSCHEME_FFF

   #include "Equations/Tools/Eigen3DTools.hpp"

   namespace GeoMHDiSCC {

      namespace Equations {

         SharedIEigenTools eigenSelector(const CouplingInformation::IndexType indexType)
         {
            if(indexType == CouplingInformation::SINGLE)
            {
               SharedEigen3DTools spTools(new Eigen3DTools());

               return spTools;

            } else
            {
               throw Exception("Unknown eigen tools type");
            }
         }

      }
   }
#endif //QUICC_SPATIALSCHEME_FFF

// Configure code to use CFT scheme
#ifdef QUICC_SPATIALSCHEME_CFT

   #include "Equations/Tools/Eigen1DTools.hpp"

   namespace GeoMHDiSCC {

      namespace Equations {

         SharedIEigenTools eigenSelector(const CouplingInformation::IndexType indexType)
         {
            if(indexType == CouplingInformation::SLOWEST_SINGLE_RHS)
            {
               SharedEigen1DTools spTools(new Eigen1DTools());

               return spTools;

            } else
            {
               throw Exception("Unknown eigen tools type");
            }
         }

      }
   }
#endif //QUICC_SPATIALSCHEME_CFT

// Configure code to use AFT scheme
#ifdef QUICC_SPATIALSCHEME_AFT

   #include "Equations/Tools/Eigen1DTools.hpp"

   namespace GeoMHDiSCC {

      namespace Equations {

         SharedIEigenTools eigenSelector(const CouplingInformation::IndexType indexType)
         {
            if(indexType == CouplingInformation::SLOWEST_SINGLE_RHS)
            {
               SharedEigen1DTools spTools(new Eigen1DTools());

               return spTools;

            } else
            {
               throw Exception("Unknown eigen tools type");
            }
         }

      }
   }
#endif //QUICC_SPATIALSCHEME_AFT

// Configure code to use BLFl, SLFl, WLFl schemes
#if defined QUICC_SPATIALSCHEME_BLFL || defined QUICC_SPATIALSCHEME_SLFL || defined QUICC_SPATIALSCHEME_WLFL

   #include "Equations/Tools/EigenSHlTools.hpp"
   #include "Equations/Tools/EigenSHlmTools.hpp"

   namespace GeoMHDiSCC {

      namespace Equations {

         SharedIEigenTools eigenSelector(const CouplingInformation::IndexType indexType)
         {
            if(indexType == CouplingInformation::SLOWEST_MULTI_RHS)
            {
               SharedEigenSHlTools spTools(new EigenSHlTools());

               return spTools;

            } else if(indexType == CouplingInformation::MODE)
            {
               SharedEigenSHlmTools spTools(new EigenSHlmTools());

               return spTools; 

            } else
            {
               throw Exception("Unknown eigen tools type");
            }
         }
      }
   }

#endif //defined QUICC_SPATIALSCHEME_BLFL || defined QUICC_SPATIALSCHEME_SLFL || defined QUICC_SPATIALSCHEME_WLFL

// Configure code to use BLFm, SLFm, WLFm schemes
#if defined QUICC_SPATIALSCHEME_BLFM || defined QUICC_SPATIALSCHEME_SLFM || defined QUICC_SPATIALSCHEME_WLFM

   #include "Equations/Tools/EigenSHmTools.hpp"

   namespace GeoMHDiSCC {

      namespace Equations {

         SharedIEigenTools eigenSelector(const CouplingInformation::IndexType indexType)
         {
            if(indexType == CouplingInformation::SLOWEST_SINGLE_RHS)
            {
               SharedEigenSHmTools spTools(new EigenSHmTools());

               return spTools;

            } else
            {
               throw Exception("Unknown eigen tools type");
            }
         }
      }
   }
#endif //defined QUICC_SPATIALSCHEME_BLFM || defined QUICC_SPATIALSCHEME_SLFM || defined QUICC_SPATIALSCHEME_WLFM

// Configure code to use WFT scheme
#ifdef QUICC_SPATIALSCHEME_WFT

   #include "Equations/Tools/Eigen1DTools.hpp"

   namespace GeoMHDiSCC {

      namespace Equations {

         SharedIEigenTools eigenSelector(const CouplingInformation::IndexType indexType)
         {
            if(indexType == CouplingInformation::SLOWEST_SINGLE_RHS)
            {
               SharedEigen1DTools spTools(new Eigen1DTools());

               return spTools;

            } else
            {
               throw Exception("Unknown eigen tools type");
            }
         }

      }
   }
#endif //QUICC_SPATIALSCHEME_WFT

// Configure code to use TT scheme
#ifdef QUICC_SPATIALSCHEME_TT
   
   #include "Equations/Tools/NoEigenTools.hpp"

   namespace GeoMHDiSCC {

      namespace Equations {

         SharedIEigenTools eigenSelector(const CouplingInformation::IndexType indexType)
         {
            if(indexType == CouplingInformation::SINGLE)
            {
               SharedNoEigenTools spTools(new NoEigenTools());

               return spTools;

            } else
            {
               throw Exception("Unknown eigen tools type");
            }
         }

      }
   }
#endif //QUICC_SPATIALSCHEME_TT

// Configure code to use TF scheme
#ifdef QUICC_SPATIALSCHEME_TF
   
   #include "Equations/Tools/Eigen1DTools.hpp"

   namespace GeoMHDiSCC {

      namespace Equations {

         SharedIEigenTools eigenSelector(const CouplingInformation::IndexType indexType)
         {
            if(indexType == CouplingInformation::SLOWEST_SINGLE_RHS)
            {
               SharedEigen1DTools spTools(new Eigen1DTools());

               return spTools;

            } else
            {
               throw Exception("Unknown eigen tools type");
            }
         }

      }
   }
#endif //QUICC_SPATIALSCHEME_TF
