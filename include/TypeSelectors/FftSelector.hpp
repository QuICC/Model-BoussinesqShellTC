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
   #include "FastTransforms/CShellChebyshevFftwTransform.hpp"
   #include "FastTransforms/SShellChebyshevFftwTransform.hpp"

   namespace GeoMHDiSCC {

      namespace Transform {

         /// Typedef for FFTW's FFT tools implementation
         typedef FftwTools FftToolsType;

         /// Typedef for FFTW's FFT implementation
         typedef FftwTransform FftTransformType;

         /// Typedef for FFTW's Chebyshev FFT implementation
         typedef ChebyshevFftwTransform ChebyshevTransformType;

         /// Typedef for FFTW's cylindrical shell Chebyshev FFT implementation
         typedef CShellChebyshevFftwTransform CShellChebyshevTransformType;

         /// Typedef for FFTW's spherical shell Chebyshev FFT implementation
         typedef SShellChebyshevFftwTransform SShellChebyshevTransformType;
      }
   }
#endif //GEOMHDISCC_FFT_FFTW

#endif // FFTSELECTOR_HPP
