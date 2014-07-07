/**
 * @file PythonCoordinator.hpp
 * @brief Small coordinator for the Python interpreter to initialise and finalize with multiple wrappers
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef PYTHONCOORDINATOR_HPP
#define PYTHONCOORDINATOR_HPP

// First include
//
#include "Python/PythonHeader.hpp"

// System includes
//

// External includes
//

// Project includes
//

namespace GeoMHDiSCC {

   /**
    * @brief Small coordinator for the Python interpreter to initialise and finalize with multiple wrappers
    */
   class PythonCoordinator
   {
      public:
         /**
          * @brief Initialise the Python interpreter
          */
         static void init();

         /**
          * @brief Register Python wrapper
          */
         static void registerWrapper();

         /**
          * @brief Unregister Python wrapper
          */
         static void unregisterWrapper();

         /**
          * @brief Finalize the Python interpreter
          */
         static void finalize();
         
      protected:

      private:
         /**
          * @brief Counter for the number of active FFTW objects
          */
         static int sCounter; 

         /**
          * @brief Constructor
          */
         PythonCoordinator();

         /**
          * @brief Destructor
          */
         ~PythonCoordinator();
   };

}

#endif // PYTHONCOORDINATOR_HPP
