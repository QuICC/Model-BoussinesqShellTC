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
#include "Exceptions/Exception.hpp"

namespace QuICC {

namespace Transform {

  // Fastest FFTW plan creation
  #ifdef QUICC_FFTPLAN_FAST
     unsigned int FftwLibrary::sPlanFlag = FFTW_ESTIMATE;
  // Medium FFTW plan creation
  #elif defined QUICC_FFTPLAN_MEDIUM
     unsigned int FftwLibrary::sPlanFlag = FFTW_MEASURE;
  // Slow FFTW plan creation
  #elif defined QUICC_FFTPLAN_SLOW
     unsigned int FftwLibrary::sPlanFlag = FFTW_PATIENT;
  #endif // QUICC_FFTW_ESTIMATE

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

   void FftwLibrary::initFft()
   {
      #if defined QUICC_THREADS_PTHREADS || defined QUICC_THREADS_OPENMP
      if(FftwLibrary::sCounter == 0)
      {
         // Initialize FFTW's threads
         int error = fftw_init_threads();

         if(error == 0)
         {
            throw Exception("FFTW's threads initialization failed!");
         }
      }

      // Number of threads
      fftw_plan_with_nthreads(2);
      #endif //defined QUICC_THREADS_PTHREADS || defined QUICC_THREADS_OPENMP
   }

   void FftwLibrary::cleanupFft()
   {
      // Check if all objects have been destroyed
      if(FftwLibrary::sCounter == 0)
      {
         #if defined QUICC_THREADS_PTHREADS || defined QUICC_THREADS_OPENMP
            fftw_cleanup_threads();
         #else
            fftw_cleanup();
         #endif //defined QUICC_THREADS_PTHREADS || defined QUICC_THREADS_OPENMP
      }
   }

}
}
