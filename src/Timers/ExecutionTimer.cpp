/** \file ExecutionTimer.cpp
 *  \brief Source of the implementation of an execution timer
 */

// Configuration includes
//
#include "Framework/FrameworkMacro.h"

// System includes
//

// External includes
//

// Class include
//
#include "Timers/ExecutionTimer.hpp"

// Project includes
//
#include "IoTools/Formatter.hpp"

namespace GeoMHDiSCC {

   ExecutionTimer::ExecutionTimer(const bool autostart)
      : TimerMacro(autostart), mTimes(Array::Zero(TOTAL + 1))
   {
   }

   ExecutionTimer::~ExecutionTimer()
   {
   }

   void ExecutionTimer::update(BreakPoint point)
   {
      // Increment measured time
      this->mTimes(point) += this->time();

      // Unless TOTAL is request, add time to total execution time
      if(point != TOTAL)
      {
         this->mTimes(TOTAL) += this->time();
      }

      std::cerr << this->mTimes.transpose() << std::endl;
   }

   void ExecutionTimer::analyze(Array& min, Array& max)
   {
      // Get the "global" times from MPI code
      #ifdef GEOMHDISCC_MPI
         // Resize the storage
         min.resize(this->mTimes.size());
         max.resize(this->mTimes.size());

         // Get the max values
         MPI_Allreduce(this->mTimes.data(), max.data(), this->mTimes.size(), MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

         // Get the min values
         MPI_Allreduce(this->mTimes.data(), min.data(), this->mTimes.size(), MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
         // Get the mean values
         MPI_Allreduce(MPI_IN_PLACE, this->mTimes.data(), this->mTimes.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

         // Compute mean times per CPU
         this->mTimes /= static_cast<MHDFloat>(FrameworkMacro::nCpu());
      #endif // GEOMHDISCC_MPI
   }

   void ExecutionTimer::printInfo(std::ostream& stream)
   {
      // Analyze the data
      Array min;
      Array max;
      this->analyze(min, max);

      // Create nice looking ouput header
      IoTools::Formatter::printNewline(stream);
      IoTools::Formatter::printLine(stream, '-');
      IoTools::Formatter::printCentered(stream, "Execution time information", '*');
      IoTools::Formatter::printLine(stream, '-');

      std::stringstream oss;

      // get a nice base for info
      int base = 20;
      oss << std::fixed << std::setprecision(1) << this->mTimes.maxCoeff();
      base += oss.str().size() + 1;
      oss.str("");

      // Output initialisation time
      oss << "Initialisation [s]: " << std::fixed << std::setprecision(1) << this->mTimes(INIT);
      if(max.size() != 0)
      {
         oss << " / " << max(INIT) << " / " << min(INIT);
      }
      IoTools::Formatter::printCentered(stream, oss.str(), ' ', base);
      oss.str("");

      // Output prerun time
      oss << "PreRun [s]: " << std::fixed << std::setprecision(1) << this->mTimes(PRERUN);
      if(max.size() != 0)
      {
         oss << " / " << max(PRERUN) << " / " << min(PRERUN);
      }
      IoTools::Formatter::printCentered(stream, oss.str(),' ', base);
      oss.str("");

      // Output computation time
      oss << "Computation [s]: " << std::fixed << std::setprecision(1) << this->mTimes(RUN);
      if(max.size() != 0)
      {
         oss << " / " << max(RUN) << " / " << min(RUN);
      }
      IoTools::Formatter::printCentered(stream, oss.str(), ' ', base);
      oss.str("");

      // Output postrun time
      oss << "PostRun [s]: " << std::fixed << std::setprecision(1) << this->mTimes(POSTRUN);
      if(max.size() != 0)
      {
         oss << " / " << max(POSTRUN) << " / " << min(POSTRUN);
      }
      IoTools::Formatter::printCentered(stream, oss.str(), ' ', base);
      oss.str("");

      IoTools::Formatter::printLine(stream, '-');

      // Output total execution time
      oss << "Total execution time: " << std::fixed << std::setprecision(1) << this->mTimes(TOTAL);
      if(max.size() != 0)
      {
         oss << " / " << max(TOTAL) << " / " << min(TOTAL);
      }
      oss << " seconds";
      IoTools::Formatter::printCentered(stream, oss.str(), ' ', base);
      oss.str("");

      IoTools::Formatter::printLine(stream, '*');
   }

}
