/**
 * @file CuFftTools.hpp
 * @brief Definition of some useful constants and tools for cuFFT 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef CUFFTTOOLS_HPP
#define CUFFTTOOLS_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"

namespace GeoMHDiSCC {

namespace Transform {

   /**
    * @brief Contains some useful constants and tools for FFTW
    */
   class CuFftTools
   {
      public:
         /**
          * @brief Compute the dealiased size for FFT (R<->R, or Z<->Z)
          *
          * @param size Size to dealias
          */
         static int dealiasFft(const int size);

         /**
          * @brief Compute the dealiased size for FFT (R<->Z, or Z<->R)
          *
          * @param size Size to dealias
          */
         static int dealiasMixedFft(const int size);

         /**
          * @brief Compute the dealiased size for cosine FFT
          *
          * @param size Size to dealias
          */
         static int dealiasCosFft(const int size);

         /**
          * @brief Optimise the FFT sizes
          *
          * @param size Current size of the FFT
          */
         static int optimizeFft(const int size);
         
      protected:

      private:
         /**
          * @brief Standard dealiasing factor (usually 3/2)
          */
         static const MHDFloat STD_DEALIASING;

         /**
          * @brief The real <-> complex fast Fourrier transform dealiasing factor
          */
         static const MHDFloat MIXED_DEALIASING;

         /**
          * @brief Cosine dealiasing factor (usually 3/2)
          */
         static const MHDFloat COS_DEALIASING;

         /**
          * @brief Maximul extension width to consider for optimization
          */
         static const MHDFloat OPTIMIZATION_WIDTH;

         /**
          * @brief Empty constructor
          */
         CuFftTools();

         /**
          * @brief Empty Destructor
          */
         ~CuFftTools();

   };

}
}

#endif // CUFFTTOOLS_HPP
