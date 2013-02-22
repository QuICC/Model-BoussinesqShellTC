/** \file FftwTools.cpp
 *  \brief Defines some useful constants and tools for FFTW
 */

// System includes
//
#include <cmath>

// External includes
//

// Class include
//
#include "FastTransforms/FftwTools.hpp"

// Project includes
//
#include "IoTools/Formatter.hpp"

namespace GeoMHDiSCC {

namespace Transform {

   const MHDFloat FftwTools::STD_DEALIASING = 3.0/2.0;

   const MHDFloat FftwTools::MIXED_DEALIASING = 3.0;

   const MHDFloat FftwTools::OPTIMIZATION_WIDTH = 0.05;

   int FftwTools::dealiasFft(const int size)
   {
      return std::ceil(FftwTools::STD_DEALIASING*static_cast<MHDFloat>(size));
   }

   int FftwTools::dealiasMixedFft(const int size)
   {
      return std::ceil(FftwTools::MIXED_DEALIASING*static_cast<MHDFloat>(size));
   }

   int FftwTools::optimiseFft(const int size)
   {
      // Store the input size
      int factorised;

      // optimisation width (maximum of 5% extension)
      int width = std::max(3, static_cast<int>(FftwTools::OPTIMIZATION_WIDTH*size));

      bool good = false;
      int opt = -1;

      // Loop while optimisation is possible
      while(!good && opt < width)
      {
         // Increment optimisation
         opt++;

         // Get trial size
         factorised = size + opt;

         // Compute factorisation (trial division algorithm)
         while(factorised > 1)
         {
            // Check for optimised factor 2
            if(factorised % 2 == 0)
            {
               factorised /= 2;
               good = true;
            // Check for optimised factor 3
            } else if(factorised % 3 == 0)
            {
               factorised /= 3;
               good = true;
            // Check for optimised factor 5
            } else if(factorised % 5 == 0)
            {
               factorised /= 5;
               good = true;
            // Check for optimised factor 7
            } else if(factorised % 7 == 0)
            {
               factorised /= 7;
               good = true;
            // Check for optimised factor 11 (single factor)
            } else if(factorised % 11 == 0)
            {
               factorised /= 11;

               // Only single factor possible
               if(factorised == 1)
               {
                  good = true;
               }
               break;
            // Check for optimised factor 13 (single factor)
            } else if(factorised % 13 == 0)
            {
               factorised /= 13;

               // Only single factor possible
               if(factorised == 1)
               {
                  good = true;
               }
               break;
            } else
            {
               good = false;
               break;
            }
         }
      }

      // Size is not good and optimisation failed
      if(!good)
      {
         // Create nice looking warnig message
         IoTools::Formatter::printLine(std::cout, '%');
         IoTools::Formatter::printCentered(std::cout, "WARNING: FFT size optimization failed", '%');
         IoTools::Formatter::printCentered(std::cout, "Selected FFT size might be slow!", '%');
         IoTools::Formatter::printLine(std::cout, '%');
         IoTools::Formatter::printNewline(std::cout);

         return size;

      // Size was not good but optimisation worked
      } else if(good && opt != 0)
      {
         // Create nice looking warning message
         IoTools::Formatter::printLine(std::cout, '%');
         std::stringstream oss;
         if(opt > 2)
         {
            oss << "WARNING: ";
         }
         oss << "Extended FFT size (+" << opt << ")!";
         IoTools::Formatter::printCentered(std::cout, oss.str(), '%');
         IoTools::Formatter::printLine(std::cout, '%');
         IoTools::Formatter::printNewline(std::cout);

         return size + opt;
      }
   }

}
}
