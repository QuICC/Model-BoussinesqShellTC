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

   Array ProfilerBase::timings = Array::Zero(static_cast<int>(ProfilerBase::NBREAKPOINT));

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
         case BWDND:
            return "Bwd ND";
         case FWD1D:
            return "Fwd 1D";
         case FWD2D:
            return "Fwd 2D";
         case FWDND:
            return "Fwd ND";

         // Below this line are the "detailed" break points level 2
         case BWDDEALIAS:
            return "Bwd dealias";
         case BWD1DTRA:
            return "Bwd 1D Transform";
         case BWD2DTRA:
            return "Bwd 2D Transform";
         case BWDNDTRA:
            return "Bwd ND Transform";
         case FWD1DTRA:
            return "Fwd 1D Transform";
         case FWD2DTRA:
            return "Fwd 2D Transform";
         case FWDNDTRA:
            return "Fwd ND Transform";
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

         // Below this line are the "detailed" break points level 3
         case BWD1DTRAFFT:
            return "Bwd 1D FFT";
         case BWD1DTRAR:
            return "Bwd 1D Solve R x = y";
         case BWD1DTRAR2:
            return "Bwd 1D Solve R^2 x = y";
         case BWD1DTRADIFF:
            return "Bwd 1D Solve D^-1 x = y";
         case BWD1DTRADIFF2:
            return "Bwd 1D Solve D^-2 x = y";
         case FWD1DTRAFFT:
            return "Fwd 1D FFT";
         case TSTEPMPI:
            return "timestep MPI";

         // Below this line are the unspecific probing break points
         case PROBEA:
            return "probe A";
         case PROBEB:
            return "probe B";
         case PROBEC:
            return "probe C";

         // Default output for unknown break point
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
