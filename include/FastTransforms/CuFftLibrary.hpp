/**
 * @file CuFftLibrary.hpp
 * @brief Static interface to the global features of the FFTW library  
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef CUFFTLIBRARY_HPP
#define CUFFTLIBRARY_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//

namespace GeoMHDiSCC {

namespace Transform {

   /**
    * @brief Static interface to the global features of the FFTW library 
    */ 
   class CuFftLibrary
   {
      public:
         /**
          * @brief Register object using FFTW library
          */
         static void registerFft();

         /**
          * @brief Unregister object using FFTW library
          */
         static void unregisterFft();

         /**
          * @brief Cleanup the FFTW library
          */
         static void cleanupFft();

      private:
         /**
          * @brief Counter for the number of active FFTW objects
          */
         static int sCounter; 

         /**
          * @brief Empty constructor
          */
         CuFftLibrary();

         /**
          * @brief Empty destructor
          */
         ~CuFftLibrary(); 
   };

}
}

#endif // CUFFTLIBRARY_HPP
