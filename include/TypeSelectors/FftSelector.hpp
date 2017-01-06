/** 
 * @file FftSelector.hpp
 * @brief Selector for the FFT implementation to be used in the code
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef FFTSELECTOR_HPP
#define FFTSELECTOR_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//

/// Setup typedef for FFTW based transforms
#if defined QUICC_FFT_FFTW

   // FFTW includes
   #include "FastTransforms/FftwTools.hpp"
   #include "FastTransforms/FftwTransform.hpp"
   #include "FastTransforms/ChebyshevFftwTransform.hpp"
   #if defined QUICC_SPATIALSCHEME_CFT
      #include "FastTransforms/CylinderChebyshevFftwTransform.hpp"
   #elif defined QUICC_SPATIALSCHEME_AFT
      #include "FastTransforms/AnnulusChebyshevFftwTransform.hpp"
   #elif defined QUICC_SPATIALSCHEME_BLFL || defined QUICC_SPATIALSCHEME_BLFM
      #include "FastTransforms/SphereChebyshevFftwTransform.hpp"
   #elif defined QUICC_SPATIALSCHEME_SLFL || defined QUICC_SPATIALSCHEME_SLFM
      #include "FastTransforms/ShellChebyshevFftwTransform.hpp"
   #endif //defined QUICC_SPATIALSCHEME_CFT

   namespace GeoMHDiSCC {

      namespace Transform {

         namespace Fft {

            /// Typedef for FFTW's FFT tools implementation
            typedef FftwTools ToolsSelector;

            /// Typedef for FFTW's FFT implementation
            typedef FftwTransform FftSelector;

            /// Typedef for FFTW's Chebyshev FFT implementation
            typedef ChebyshevFftwTransform ChebyshevSelector;

            #if defined QUICC_SPATIALSCHEME_CFT
               /// Typedef for FFTW's cylinder Chebyshev FFT implementation
               typedef CylinderChebyshevFftwTransform CylinderChebyshevSelector;

            #elif defined QUICC_SPATIALSCHEME_AFT
               /// Typedef for FFTW's annulus Chebyshev FFT implementation
               typedef AnnulusChebyshevFftwTransform AnnulusChebyshevSelector;

            #elif defined QUICC_SPATIALSCHEME_BLFL || defined QUICC_SPATIALSCHEME_BLFM
               /// Typedef for FFTW's sphere Chebyshev FFT implementation
               typedef SphereChebyshevFftwTransform SphereChebyshevSelector;

            #elif defined QUICC_SPATIALSCHEME_SLFL || defined QUICC_SPATIALSCHEME_SLFM
               /// Typedef for FFTW's shell Chebyshev FFT implementation
               typedef ShellChebyshevFftwTransform ShellChebyshevSelector;
            #endif //defined QUICC_SPATIALSCHEME_CFT
         }
      }
   }
#elif defined QUICC_FFT_CUFFT

   // FFTW includes
   #include "FastTransforms/CuFftTools.hpp"
   #include "FastTransforms/CuFftTransform.hpp"
   #include "FastTransforms/ChebyshevCuFftTransform.hpp"

   namespace GeoMHDiSCC {

      namespace Transform {

         namespace Fft {

            /// Typedef for FFTW's FFT tools implementation
            typedef CuFftTools ToolsSelector;

            /// Typedef for FFTW's FFT implementation
            typedef CuFftTransform FftSelector;

            /// Typedef for FFTW's Chebyshev FFT implementation
            typedef ChebyshevCuFftTransform ChebyshevSelector;
         }
      }
   }

#endif //QUICC_FFT_FFTW

#endif // FFTSELECTOR_HPP
