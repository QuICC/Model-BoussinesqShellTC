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
         enum BreakPoint {BWDTRANSFORM, NONLINEAR, FWDTRANSFORM, PROGNOSTICEQUATION, DIAGNOSTICEQUATION, TRIVIALEQUATION, CONTROL, IO, BWD1D, BWD2D, BWD3D, FWD1D, FWD2D, FWD3D, FWDCONVSEND, BWDCONVSEND, FWDCONVRECV, BWDCONVRECV, TSTEPIN, TSTEPRHS, TSTEPSOLVE, TSTEPOUT};

         /**
          * @brief Stores the maximum break point ID depending on configuration
          */
         static const int NBREAKPOINT;

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
