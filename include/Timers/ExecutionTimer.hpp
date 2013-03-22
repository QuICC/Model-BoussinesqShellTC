/** \file ExecutionTimer.hpp
 *  \brief Implementation of an execution timer
 */

#ifndef EXECUTIONTIMER_HPP
#define EXECUTIONTIMER_HPP

// Configuration includes
//
#include "Timers/TimerMacro.h"

// System includes
//

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"

namespace GeoMHDiSCC {

   /**
    * \brief Implementation of an execution timer
    */
   class ExecutionTimer: public TimerMacro
   {
      public:
         /**
          * @brief Create list of possible timing intervals
          */
         enum BreakPoint {INIT, PRERUN, RUN, POSTRUN, TOTAL};

         /**
          * @brief Constructor
          *
          * @param autostart Should the timer start at creation ?
          */
         explicit ExecutionTimer(const bool autostart = false);

         /**
          * @brief Destructor
          */
         ~ExecutionTimer();

         /**
          * @brief Update the timing of given region. 
          *
          * The possible breakpoints are defined in the BreakPoint enum
          *
          * @param point Breakpoint for which the timing has to be updated
          */
         void update(BreakPoint point);
         
         /**
          * @brief Print execution time information to stream
          *
          * @param stream  Output stream
          */
         void printInfo(std::ostream& stream);

      protected:

      private:
         /**
          * @brief Storage for the execution times
          */
         Array mTimes;

         /**
          * @brief Analyze the measured times
          */
         void analyze(Array& min, Array& max);
   };

}

#endif // EXECUTIONTIMER_HPP
