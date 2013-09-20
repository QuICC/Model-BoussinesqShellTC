/** 
 * @file SpatialSchemeTestEquationSelector.hpp
 * @brief Implementation of test equation selector for a spatial scheme
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef SPATIALSCHEMETESTEQUATIONSELECTOR_HPP
#define SPATIALSCHEMETESTEQUATIONSELECTOR_HPP

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
   #include "Equations/Tests/TestTTTForwardScalar.hpp"
   #include "Equations/Tests/TestTTTBackwardScalar.hpp"

   namespace GeoMHDiSCC {

      namespace TestSuite {

         /**
          * @brief Specialisation for TTT scheme
          */
         struct SpatialSchemeTestEquationSelector
         {
            /**
             * @brief Forward scalar test
             */
            typedef Equations::TestTTTForwardScalar ForwardScalarType;

            /**
             * @brief Backward scalar test
             */
            typedef Equations::TestTTTBackwardScalar BackwardScalarType;
         };

      }
   }
#endif //GEOMHDISCC_SPATIALSCHEME_TTT

// Configure code to use TFT scheme
#ifdef GEOMHDISCC_SPATIALSCHEME_TFT

   // TFT includes
   #include "Equations/Tests/TestTFTForwardScalar.hpp"
   #include "Equations/Tests/TestTFTBackwardScalar.hpp"

   namespace GeoMHDiSCC {

      namespace TestSuite {

         /**
          * @brief Specialisation for TFT scheme
          */
         struct SpatialSchemeTestEquationSelector
         {
            /**
             * @brief Forward scalar test
             */
            typedef Equations::TestTFTForwardScalar ForwardScalarType;

            /**
             * @brief Backward scalar test
             */
            typedef Equations::TestTFTBackwardScalar BackwardScalarType;
         };

      }
   }
#endif //GEOMHDISCC_SPATIALSCHEME_TFT

// Configure code to use TFF scheme
#ifdef GEOMHDISCC_SPATIALSCHEME_TFF

   // TFF includes
   #include "Equations/Tests/TestTFFForwardScalar.hpp"
   #include "Equations/Tests/TestTFFBackwardScalar.hpp"

   namespace GeoMHDiSCC {

      namespace TestSuite {

         /**
          * @brief Specialisation for TFF scheme
          */
         struct SpatialSchemeTestEquationSelector
         {
            /**
             * @brief Forward scalar test
             */
            typedef Equations::TestTFFForwardScalar ForwardScalarType;

            /**
             * @brief Backward scalar test
             */
            typedef Equations::TestTFFBackwardScalar BackwardScalarType;
         };

      }
   }
#endif //GEOMHDISCC_SPATIALSCHEME_TFF

// Configure code to use FFF scheme
#ifdef GEOMHDISCC_SPATIALSCHEME_FFF

   // FFF includes
   #include "Equations/Tests/TestFFFForwardScalar.hpp"
   #include "Equations/Tests/TestFFFBackwardScalar.hpp"

   namespace GeoMHDiSCC {

      namespace TestSuite {

         /**
          * @brief Specialisation for FFF scheme
          */
         struct SpatialSchemeTestEquationSelector
         {
            /**
             * @brief Forward scalar test
             */
            typedef Equations::TestFFFForwardScalar ForwardScalarType;

            /**
             * @brief Backward scalar test
             */
            typedef Equations::TestFFFBackwardScalar BackwardScalarType;
         };

      }
   }
#endif //GEOMHDISCC_SPATIALSCHEME_FFF

// Configure code to use CFT scheme
#ifdef GEOMHDISCC_SPATIALSCHEME_CFT

   // CFT includes
   #include "Equations/Tests/TestCFTForwardScalar.hpp"
   #include "Equations/Tests/TestCFTBackwardScalar.hpp"

   namespace GeoMHDiSCC {

      namespace TestSuite {

         /**
          * @brief Specialisation for CFT scheme
          */
         struct SpatialSchemeTestEquationSelector
         {
            /**
             * @brief Forward scalar test
             */
            typedef Equations::TestCFTForwardScalar ForwardScalarType;

            /**
             * @brief Backward scalar test
             */
            typedef Equations::TestCFTBackwardScalar BackwardScalarType;
         };

      }
   }
#endif //GEOMHDISCC_SPATIALSCHEME_CFT

// Configure code to use SLF scheme
#ifdef GEOMHDISCC_SPATIALSCHEME_SLF

   // SLF includes
   #include "Equations/Tests/TestSLFForwardScalar.hpp"
   #include "Equations/Tests/TestSLFBackwardScalar.hpp"

   namespace GeoMHDiSCC {

      namespace TestSuite {

         /**
          * @brief Specialisation for SLF scheme
          */
         struct SpatialSchemeTestEquationSelector
         {
            /**
             * @brief Forward scalar test
             */
            typedef Equations::TestSLFForwardScalar ForwardScalarType;

            /**
             * @brief Backward scalar test
             */
            typedef Equations::TestSLFBackwardScalar BackwardScalarType;
         };

      }
   }
#endif //GEOMHDISCC_SPATIALSCHEME_SLF

// Configure code to use WFT scheme
#ifdef GEOMHDISCC_SPATIALSCHEME_WFT

   // WFT includes
   #include "Equations/Tests/TestWFTForwardScalar.hpp"
   #include "Equations/Tests/TestWFTBackwardScalar.hpp"

   namespace GeoMHDiSCC {

      namespace TestSuite {

         /**
          * @brief Specialisation for WFT scheme
          */
         struct SpatialSchemeTestEquationSelector
         {
            /**
             * @brief Forward scalar test
             */
            typedef Equations::TestWFTForwardScalar ForwardScalarType;

            /**
             * @brief Backward scalar test
             */
            typedef Equations::TestWFTBackwardScalar BackwardScalarType;
         };

      }
   }
#endif //GEOMHDISCC_SPATIALSCHEME_WFT

// Configure code to use WLF scheme
#ifdef GEOMHDISCC_SPATIALSCHEME_WLF

   // WLF includes
   #include "Equations/Tests/TestWLFForwardScalar.hpp"
   #include "Equations/Tests/TestWLFBackwardScalar.hpp"

   namespace GeoMHDiSCC {

      namespace TestSuite {

         /**
          * @brief Specialisation for WLF scheme
          */
         struct SpatialSchemeTestEquationSelector
         {
            /**
             * @brief Forward scalar test
             */
            typedef Equations::TestWLFForwardScalar ForwardScalarType;

            /**
             * @brief Backward scalar test
             */
            typedef Equations::TestWLFBackwardScalar BackwardScalarType;
         };

      }
   }
#endif //GEOMHDISCC_SPATIALSCHEME_WLF

#endif // SPATIALSCHEMETESTEQUATIONSELECTOR_HPP
