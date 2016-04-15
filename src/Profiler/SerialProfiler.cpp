/** 
 * @file SerialProfiler.cpp
 * @brief Source of the serial profiler implementation
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// System includes
//

// External includes
//

// Class include
//
#include "Profiler/SerialProfiler.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Debug {

   std::map<ProfilerBase::BreakPoint, timespec> SerialProfiler::t_starts = std::map<ProfilerBase::BreakPoint, timespec>();

   std::map<ProfilerBase::BreakPoint, timespec> SerialProfiler::t_stops = std::map<ProfilerBase::BreakPoint, timespec>();

   void SerialProfiler::init()
   {
      for(int i = 0; i < static_cast<int>(ProfilerBase::NBREAKPOINT); i++)
      {
         t_starts.insert(std::make_pair(static_cast<ProfilerBase::BreakPoint>(i), timespec()));
         t_stops.insert(std::make_pair(static_cast<ProfilerBase::BreakPoint>(i), timespec()));
      }
   }

   void SerialProfiler::reset()
   {
      ProfilerBase::reset();

      for(std::map<ProfilerBase::BreakPoint, timespec>::iterator it = t_starts.begin(); it != t_starts.end(); ++it)
      {
         it->second = timespec();
      }

      for(std::map<ProfilerBase::BreakPoint, timespec>::iterator it = t_stops.begin(); it != t_stops.end(); ++it)
      {
         it->second = timespec();
      }
   }

   void SerialProfiler::start(ProfilerBase::BreakPoint point)
   {
      // Get the starting timespec
      clock_gettime(CLOCK_REALTIME, &t_starts.at(point));
   }

   void SerialProfiler::stop(ProfilerBase::BreakPoint point)
   {
      // Get the stopping timespec
      clock_gettime(CLOCK_REALTIME, &t_stops.at(point));

      // Store time
      ProfilerBase::update(point, SerialProfiler::elapsedSeconds(t_starts.at(point), t_stops.at(point)));
   }

   MHDFloat SerialProfiler::elapsedSeconds(timespec &t1, timespec &t2)
   {
      // Compute elapsed seconds between the two timespecs
      return static_cast<MHDFloat>(t2.tv_sec - t1.tv_sec) + static_cast<MHDFloat>(t2.tv_nsec - t1.tv_nsec)/1.0e9;
   }

   void SerialProfiler::getTimings(Array& ts)
   {
      ts.resize(static_cast<int>(ProfilerBase::NBREAKPOINT));
      for(int i = 0; i < ts.size(); ++i)
      {
         ts(i) = ProfilerBase::time(static_cast<ProfilerBase::BreakPoint>(i));
      }
   }

   void SerialProfiler::analyze(Array& ts, Array& min, Array& max)
   {
      SerialProfiler::getTimings(ts);

      min.resize(0);
      max.resize(0);
   }

}
}
