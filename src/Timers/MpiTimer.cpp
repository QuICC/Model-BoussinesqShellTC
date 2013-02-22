/** \file MpiTimer.cpp
 *  \brief Source of the implementation of MPI timer
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
#include "Timers/MpiTimer.hpp"

// Project includes
//

namespace GeoMHDiSCC {

   MpiTimer::MpiTimer(const bool autostart)
      : mStart(0.0), mStop(0.0)
   {
      // Check if timer should be started at creation
      if(autostart)
      {
         this->start();
      }
   }

   MpiTimer::~MpiTimer()
   {
   }

   void MpiTimer::start()
   {
      // Get starting MPI time
      this->mStart = MPI_Wtime();
   }

   void MpiTimer::stop()
   {
      // Get stoppping MPI time
      this->mStop = MPI_Wtime();
   }

   MHDFloat MpiTimer::time() const
   {
      // Compute elapsed time
      return this->mStop - this->mStart;
   }

   MHDFloat MpiTimer::resetTimer()
   {
      // Stop the timer
      this->stop();

      // Get elapsed time
      MHDFloat tmp = this->time();

      // Set start time to stopping time
      this->mStart = this->mStop;

      // return elapsed time
      return tmp;
   }

}
