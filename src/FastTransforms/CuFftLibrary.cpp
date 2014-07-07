/** 
 * @file CuFftLibrary.cpp
 * @brief Source for the static interface to cuFFT
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// System includes
//

// External includes
//
#include <cufftw.h>

// Class include
//
#include "FastTransforms/CuFftLibrary.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Transform {

  // Fastest FFTW plan creation
  #ifdef GEOMHDISCC_FFTPLAN_FAST
     unsigned int CuFftLibrary::sPlanFlag = FFTW_ESTIMATE;
  // Medium FFTW plan creation
  #elif defined GEOMHDISCC_FFTPLAN_MEDIUM
     unsigned int CuFftLibrary::sPlanFlag = FFTW_MEASURE;
  // Slow FFTW plan creation
  #elif defined GEOMHDISCC_FFTPLAN_SLOW
     unsigned int CuFftLibrary::sPlanFlag = FFTW_PATIENT;
  #endif // GEOMHDISCC_FFTW_ESTIMATE

   int CuFftLibrary::sCounter = 0;

   CuFftLibrary::CuFftLibrary()
   {
   }

   CuFftLibrary::~CuFftLibrary()
   {
   }

   unsigned int CuFftLibrary::planFlag()
   {
      return CuFftLibrary::sPlanFlag;
   }

   void CuFftLibrary::registerFft()
   {
      ++CuFftLibrary::sCounter;
   }

   void CuFftLibrary::unregisterFft()
   {
      --CuFftLibrary::sCounter;
   }

   void CuFftLibrary::cleanupFft()
   {
      // Check if all objects have been destroyed
      if(CuFftLibrary::sCounter == 0)
      {
         fftw_cleanup();
      }
   }

}
}
