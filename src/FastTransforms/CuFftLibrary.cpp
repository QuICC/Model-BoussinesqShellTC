/** 
 * @file CuFftLibrary.cpp
 * @brief Source for the static interface to cuFFT
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// System includes
//

// External includes
//
#include <cufft.h>

// Class include
//
#include "FastTransforms/CuFftLibrary.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Transform {

   int CuFftLibrary::sCounter = 0;

   CuFftLibrary::CuFftLibrary()
   {
   }

   CuFftLibrary::~CuFftLibrary()
   {
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
   }

}
}
