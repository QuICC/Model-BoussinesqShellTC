/**
 * @file ProfilerBase.hpp
 * @brief Implementation of the base of a profiling timer 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef PROFILERBASE_HPP
#define PROFILERBASE_HPP

// Configuration includes
//

// System includes
//
#include <map>

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"

namespace GeoMHDiSCC {

namespace Debug {

   /**
    * @brief Implementation of the base of a profiling timer
    */
   class ProfilerBase
   {
      public:
         /**
          * @brief List of break points for the profiling
          */
         enum BreakPoint {
            // Coarse profiling points (level 0)
            BWDTRANSFORM = 0,
            NONLINEAR,
            FWDTRANSFORM,
            PROGNOSTICEQUATION,
            DIAGNOSTICEQUATION,
            TRIVIALEQUATION,
            CONTROL,
            IO,
            // Detailed profiling points level 1 (included in level 0)
            BWD1D,
            BWD2D,
            BWDND,
            FWD1D,
            FWD2D,
            FWDND,
            // Detailed profiling points level 2 (included in level 1)
            BWDDEALIAS,
            BWD1DTRA,
            BWD2DTRA,
            BWDNDTRA,
            FWD1DTRA,
            FWD2DTRA,
            FWDNDTRA,
            BWDSENDWAIT,
            BWDSENDCONV,
            BWDRECVWAIT,
            BWDRECVCONV,
            FWDSENDWAIT,
            FWDSENDCONV,
            FWDRECVWAIT,
            FWDRECVCONV,
            TSTEPIN,
            TSTEPRHS,
            TSTEPSOLVE,
            TSTEPOUT,
            // Detailed profiling points level 3 (included in level 2)
            BWD1DTRAFFT,
            BWD1DTRAR,
            BWD1DTRAR2,
            BWD1DTRADIFF,
            BWD1DTRADIFF2,
            FWD1DTRAFFT,
            // Break point bounding value
            #ifdef GEOMHDISCC_PROFILER_DETAILED
            NBREAKPOINT
            #else
            NBREAKPOINT = IO + 1
            #endif // GEOMHDISCC_PROFILER_DETAILED
         };

         /**
          * @brief Get elapsed time for provided break point
          *
          * @param point Break point
          */
         static MHDFloat time(BreakPoint point);

         /**
          * @brief Reset the profiling timings
          */
         static void reset();

         /**
          * @brief Get a human readable name for break point
          *
          * @param point Break point
          */
         static std::string pointName(BreakPoint point);
         
      protected:
         /**
          * @brief Constructor
          */
         ProfilerBase();

         /**
          * @brief Destructor
          */
         virtual ~ProfilerBase();

         /**
          * @brief Update measured time
          *
          * @param point   Location that has been profiled
          * @param time    The measured time
          */
         static void update(BreakPoint point, MHDFloat time);

      private:
         /**
          * @brief Profiling timings
          */
         static Array timings;
   };

}
}

#endif // PROFILERBASE_HPP
