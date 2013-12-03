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
#ifdef GEOMHDISCC_FFT_FFTW

   // FFTW includes
   #include "FastTransforms/FftwTools.hpp"
   #include "FastTransforms/FftwTransform.hpp"
   #include "FastTransforms/ChebyshevFftwTransform.hpp"
   #include "FastTransforms/AnnulusChebyshevFftwTransform.hpp"
   #include "FastTransforms/ShellChebyshevFftwTransform.hpp"

   namespace GeoMHDiSCC {

      namespace Transform {

         namespace Fft {

            /// Typedef for FFTW's FFT tools implementation
            typedef FftwTools ToolsSelector;

            /// Typedef for FFTW's FFT implementation
            typedef FftwTransform FftSelector;

            /// Typedef for FFTW's Chebyshev FFT implementation
            typedef ChebyshevFftwTransform ChebyshevSelector;

            /// Typedef for FFTW's cylindrical shell Chebyshev FFT implementation
            typedef AnnulusChebyshevFftwTransform AnnulusChebyshevSelector;

            /// Typedef for FFTW's spherical shell Chebyshev FFT implementation
            typedef ShellChebyshevFftwTransform ShellChebyshevSelector;
         }
      }
   }
#endif //GEOMHDISCC_FFT_FFTW

#endif // FFTSELECTOR_HPP
