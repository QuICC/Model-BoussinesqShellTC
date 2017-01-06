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

namespace QuICC {

namespace Transform {

   const int CuFftLibrary::NSTREAMS = 10;

   int CuFftLibrary::sCounter = 0;

   std::vector<cudaStream_t>  CuFftLibrary::sStream = std::vector<cudaStream_t>();

   CuFftLibrary::CuFftLibrary()
   {
   }

   CuFftLibrary::~CuFftLibrary()
   {
   }

   void CuFftLibrary::registerFft()
   {
      ++CuFftLibrary::sCounter;

      if(CuFftLibrary::sCounter == 1)
      {
         for(int i = 0; i < CuFftLibrary::NSTREAMS; ++i)
         {
            CuFftLibrary::sStream.push_back(cudaStream_t());
            checkCudaErrors(cudaStreamCreate(&CuFftLibrary::sStream.at(i)));
         }
      }
   }

   void CuFftLibrary::unregisterFft()
   {
      --CuFftLibrary::sCounter;
   }

   void CuFftLibrary::cleanupFft()
   {
      if(CuFftLibrary::sCounter == 0)
      {
         for(int i = 0; i < CuFftLibrary::NSTREAMS; ++i)
         {
            checkCudaErrors(cudaStreamDestroy(CuFftLibrary::sStream.at(i)));
         }
      }
   }

}
}
