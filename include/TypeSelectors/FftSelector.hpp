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
   #include "FastTransforms/CylindricalChebyshevFftwTransform.hpp"
   #include "FastTransforms/SphericalChebyshevFftwTransform.hpp"

   namespace GeoMHDiSCC {

      namespace Transform {

         /// Typedef for FFTW's FFT tools implementation
         typedef FftwTools FftToolsType;

         /// Typedef for FFTW's FFT implementation
         typedef FftwTransform FftTransformType;

         /// Typedef for FFTW's Chebyshev FFT implementation
         typedef ChebyshevFftwTransform ChebyshevTransformType;

         /// Typedef for FFTW's cylindrical Chebyshev FFT implementation
         typedef CylindricalChebyshevFftwTransform CylindricalChebyshevTransformType;

         /// Typedef for FFTW's spherical Chebyshev FFT implementation
         typedef SphericalChebyshevFftwTransform SphericalChebyshevTransformType;
      }
   }
#endif //GEOMHDISCC_FFT_FFTW

#endif // FFTSELECTOR_HPP
