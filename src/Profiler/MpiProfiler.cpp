/** 
 * @file MpiProfiler.cpp
 * @brief Source of the implementation of MPI profiler
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
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
#include "Profiler/MpiProfiler.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Debug {

   std::map<ProfilerBase::BreakPoint, MHDFloat> MpiProfiler::t_starts = std::map<ProfilerBase::BreakPoint, MHDFloat>();

   std::map<ProfilerBase::BreakPoint, MHDFloat> MpiProfiler::t_stops = std::map<ProfilerBase::BreakPoint, MHDFloat>();

   void MpiProfiler::init()
   {
      for(int i = 0; i < static_cast<int>(ProfilerBase::NBREAKPOINT); i++)
      {
         t_starts.insert(std::make_pair(static_cast<ProfilerBase::BreakPoint>(i), 0));
         t_stops.insert(std::make_pair(static_cast<ProfilerBase::BreakPoint>(i), 0));
      }
   }

   void MpiProfiler::reset()
   {
      ProfilerBase::reset();

      for(std::map<ProfilerBase::BreakPoint, MHDFloat>::iterator it = t_starts.begin(); it != t_starts.end(); ++it)
      {
         it->second = 0;
      }

      for(std::map<ProfilerBase::BreakPoint, MHDFloat>::iterator it = t_stops.begin(); it != t_stops.end(); ++it)
      {
         it->second = 0;
      }
   }

   void MpiProfiler::start(ProfilerBase::BreakPoint point)
   {
      // Get starting MPI time
      t_starts.at(point) = MPI_Wtime();
   }

   void MpiProfiler::stop(ProfilerBase::BreakPoint point)
   {
      // Get stoppping MPI time
      t_stops.at(point) = MPI_Wtime();

      // Store time
      ProfilerBase::update(point, t_stops.at(point)-t_starts.at(point));
   }

   void MpiProfiler::getTimings(Array& ts)
   {
      ts.resize(static_cast<int>(ProfilerBase::NBREAKPOINT));
      for(int i = 0; i < ts.size(); ++i)
      {
         ts(i) = ProfilerBase::time(static_cast<ProfilerBase::BreakPoint>(i));
      }
   }

   void MpiProfiler::analyze(Array& ts, Array& min, Array& max)
   {
      MpiProfiler::getTimings(ts);

      // Resize the storage
      min.resize(ts.size());
      max.resize(ts.size());

      // Get the max values
      MPI_Allreduce(ts.data(), max.data(), ts.size(), MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

      // Get the min values
      MPI_Allreduce(ts.data(), min.data(), ts.size(), MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

      // Get the mean values
      MPI_Allreduce(MPI_IN_PLACE, ts.data(), ts.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      int size;
      // Get the number of CPUs involved
      MPI_Comm_size(MPI_COMM_WORLD, &size);
      // Compute mean times per CPU
      ts /= static_cast<MHDFloat>(size);
   }

}
}
