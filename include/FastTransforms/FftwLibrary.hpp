/**
 * @file FftwLibrary.hpp
 * @brief Static interface to the global features of the FFTW library  
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef FFTWLIBRARY_HPP
#define FFTWLIBRARY_HPP

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
   class FftwLibrary
   {
      public:
         /**
          * @brief Get the plan flag
          */
         static unsigned int planFlag();

         /**
          * @brief Initialize the FFTW library
          */
         static void initFft();

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
          * @brief FFTW3 flags for the plan setup
          */
         static unsigned int  sPlanFlag;

         /**
          * @brief Empty constructor
          */
         FftwLibrary();

         /**
          * @brief Empty destructor
          */
         ~FftwLibrary(); 
   };

}
}

#endif // FFTWLIBRARY_HPP
