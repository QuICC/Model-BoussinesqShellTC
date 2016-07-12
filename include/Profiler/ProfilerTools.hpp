/**
 * @file ProfilerTools.hpp
 * @brief Implementation of some tools for the profiling timer 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef PROFILERTOOLS_HPP
#define PROFILERTOOLS_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//

namespace GeoMHDiSCC {

namespace Debug {

   /**
    * @brief Implementation of some tools for the profiling timer
    */
   class ProfilerTools
   {
      public:
         /**
          * @brief Write profiling output
          */
         static void writeTimings();
         
         /**
          * @brief Print profiling output
          */
         static void printInfo();
         
      protected:

      private:
         /**
          * @brief Constructor
          */
         ProfilerTools();

         /**
          * @brief Destructor
          */
         ~ProfilerTools();
   };

}
}

#endif // PROFILERTOOLS_HPP
