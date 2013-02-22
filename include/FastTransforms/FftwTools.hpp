/** \file FftwTools.hpp
 *  \brief Definition of some useful constants and tools for FFTW
 */

#ifndef FFTWTOOLS_HPP
#define FFTWTOOLS_HPP

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
   class FftwTools
   {
      public:
         /***
          * @brief Compute the dealiased size for FFT (R<->R, or Z<->Z)
          *
          * @param size Size to dealias
          */
         static int dealiasFft()

         /***
          * @brief Compute the dealiased size for FFT (R<->Z, or Z<->R)
          *
          * @param size Size to dealias
          */
         static int dealiasMixedFft(const int size)

         /**
          * @brief Optimise the FFT sizes
          *
          * @param size On input current size/on output optimised size
          */
         static int optimiseFft(const int size);
         
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
          * @brief Maximul extension width to consider for optimization
          */
         static const MHDFloat OPTIMIZATION_WIDTH;

         /**
         * @brief Empty constructor
         */
         FftwTools() {};

         /**
         * @brief Empty Destructor
         */
         ~FftwTools() {};

   };

}
}

#endif // FFTWTOOLS_HPP
