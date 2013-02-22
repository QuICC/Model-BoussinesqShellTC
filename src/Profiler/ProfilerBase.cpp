/** \file ProfilerBase.cpp
 *  \brief Source of the implementation of a profiling timer
 */

// Configuration includes
//

// System includes
//
#include <mpi.h>

// External includes
//

// Class include
//
#include "Profiler/ProfilerBase.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Debug {

   #ifdef GEOMHDISCC_PROFILER_DETAILED
      const int ProfilerBase::NBREAKPOINT = ProfilerBase::BWDCONVRECV;
   #else
      const int ProfilerBase::NBREAKPOINT = ProfilerBase::IO + 1;
   #endif // GEOMHDISCC_PROFILER_DETAILED

   Array ProfilerBase::timings = Array::Zero(ProfilerBase::NBREAKPOINT);

   MHDFloat ProfilerBase::time(ProfilerBase::BreakPoint point)
   {
      return ProfilerBase::timings(point);
   }

   void ProfilerBase::update(ProfilerBase::BreakPoint point, MHDFloat time)
   {
      // Increment measured time
      ProfilerBase::timings(point) += time;
   }

   std::string ProfilerBase::pointName(ProfilerBase::BreakPoint point)
   {
      switch(point)
      {
         case BWDTRANSFORM:
            return "Backward";
         case NONLINEAR:
            return "Nonlinear";
         case FWDTRANSFORM: 
            return "Forward";
         case LINEAR:
            return "Linear";
         case TIMESTEP:
            return "Timestep";
         case CONTROL:
            return "Control";
         case IO:
            return "IO";

         // Below this line are the "detailed" break points
         case BWD1D:
            return "Bwd 1D";
         case BWD2D:
            return "Bwd 2D";
         case BWD3D:
            return "Bwd 3D";
         case FWD3D:
            return "Fwd 1D";
         case FWD2D:
            return "Fwd 2D";
         case FWD1D:
            return "Fwd 3D";
         case BWDCONVSEND: 
            return "send 2D";
         case BWDCONVRECV:
            return "receive 2D";
         case FWDCONVSEND:
            return "send 3D";
         case FWDCONVRECV:
            return "receive 3D";
      }
   }

   void ProfilerBase::reset()
   {
      ProfilerBase::timings.setZero();
   }

}
}