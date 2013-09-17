/** 
 * @file FftwLibrary.cpp
 * @brief Source for the static interface to FFTW
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// System includes
//

// External includes
//
#include <fftw3.h>

// Class include
//
#include "FastTransforms/FftwLibrary.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Transform {

  // Fastest FFTW plan creation
  #ifdef GEOMHDISCC_FFTPLAN_FAST
     unsigned int FftwLibrary::sPlanFlag = FFTW_ESTIMATE;
  // Medium FFTW plan creation
  #elif defined GEOMHDISCC_FFTPLAN_MEDIUM
     unsigned int FftwLibrary::sPlanFlag = FFTW_MEASURE;
  // Slow FFTW plan creation
  #elif defined GEOMHDISCC_FFTPLAN_SLOW
     unsigned int FftwLibrary::sPlanFlag = FFTW_PATIENT;
  #endif // GEOMHDISCC_FFTW_ESTIMATE

   int FftwLibrary::sCounter = 0;

   FftwLibrary::FftwLibrary()
   {
   }

   FftwLibrary::~FftwLibrary()
   {
   }

   unsigned int FftwLibrary::planFlag()
   {
      return FftwLibrary::sPlanFlag;
   }

   void FftwLibrary::registerFft()
   {
      ++FftwLibrary::sCounter;
   }

   void FftwLibrary::unregisterFft()
   {
      --FftwLibrary::sCounter;
   }

   void FftwLibrary::cleanupFft()
   {
      // Check if all objects have been destroyed
      if(FftwLibrary::sCounter == 0)
      {
         fftw_cleanup();
      }
   }

}
}
