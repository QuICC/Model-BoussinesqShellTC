/** \file ProfilerTools.hpp
 *  \brief Implementation of some tools for the profiling timer
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
    * \brief Implementation of the base of a profiling timer
    */
   class ProfilerTools
   {
      public:
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
