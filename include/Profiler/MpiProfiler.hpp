/** \file MpiProfiler.hpp
 *  \brief Implementation of a MPI profiling timer
 */

#ifndef MPIPROFILER_HPP
#define MPIPROFILER_HPP

// System includes
//
#include <map>

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Profiler/ProfilerBase.hpp"

namespace GeoMHDiSCC {

namespace Debug {

   /**
    * \brief Implementation of a MPI timer
    */
   class MpiProfiler: public ProfilerBase
   {
      public:
         /**
          * @brief initialise the timers
          */
         static void init();

         /**
          * @brief Start clock
          *
          * @param point Location that has been profiled
          */
         static void start(ProfilerBase::BreakPoint point);

         /**
          * @brief Stop clock
          *
          * @param point Location that has been profiled
          */
         static void stop(ProfilerBase::BreakPoint point);

         /**
          * @brief Analyze the measured times among whole framework
          *
          * @param ts   Timings for the breakpoints
          * @param min  Minimal value within framework
          * @param max  Maximal value within framework
          */
         static void analyze(Array& ts, Array& min, Array& max);
         
      protected:
         /**
          * @brief Constructor
          */
         MpiProfiler() {};

         /**
          * @brief Destructor
          */
         ~MpiProfiler() {};

      private:
         /**
          * @brief Start
          */
         static std::map<ProfilerBase::BreakPoint, MHDFloat> t_starts;

         /**
          * @brief Stop
          */
         static std::map<ProfilerBase::BreakPoint, MHDFloat> t_stops;
   };

}
}

#endif // MPIPROFILER_HPP