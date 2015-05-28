/** 
 * @file ProfilerBase.cpp
 * @brief Source of the implementation of a profiling timer
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Configuration includes
//

// System includes
//

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
      const int ProfilerBase::NBREAKPOINT = ProfilerBase::TSTEPOUT + 1;
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
         case PROGNOSTICEQUATION:
            return "Prognostic";
         case DIAGNOSTICEQUATION:
            return "Diagnostic";
         case TRIVIALEQUATION:
            return "Trivial";
         case CONTROL:
            return "Control";
         case IO:
            return "IO";

         // Below this line are the "detailed" break points level 1
         case BWD1D:
            return "Bwd 1D";
         case BWD2D:
            return "Bwd 2D";
         case BWD3D:
            return "Bwd 3D";
         case FWD1D:
            return "Fwd 1D";
         case FWD2D:
            return "Fwd 2D";
         case FWD3D:
            return "Fwd 3D";

         // Below this line are the "detailed" break points level 2
         case BWD1DTRA:
            return "Bwd 1D Transform";
         case BWD2DTRA:
            return "Bwd 2D Transform";
         case BWD3DTRA:
            return "Bwd 3D Transform";
         case FWD1DTRA:
            return "Fwd 1D Transform";
         case FWD2DTRA:
            return "Fwd 2D Transform";
         case FWD3DTRA:
            return "Fwd 3D Transform";
         case BWDSENDWAIT: 
            return "Send BWD: wait";
         case BWDSENDCONV: 
            return "Send BWD: convert";
         case BWDRECVWAIT:
            return "Receive BWD: wait";
         case BWDRECVCONV:
            return "Receive BWD: convert";
         case FWDSENDWAIT:
            return "Send FWD: wait";
         case FWDSENDCONV:
            return "Send FWD: convert";
         case FWDRECVWAIT:
            return "Receive FWD: wait";
         case FWDRECVCONV:
            return "Receive FWD: convert";
         case TSTEPIN:
            return "timestep input";
         case TSTEPRHS:
            return "timestep rhs";
         case TSTEPSOLVE:
            return "timestep solve";
         case TSTEPOUT:
            return "timestep output";
         default:
            return "Unknown break point";
      }
   }

   void ProfilerBase::reset()
   {
      ProfilerBase::timings.setZero();
   }

}
}
