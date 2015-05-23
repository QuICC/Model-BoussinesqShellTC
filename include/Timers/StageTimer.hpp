/**
 * @file StageTimer.hpp
 * @brief Implementation of a very simple stage timer 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef STAGETIMER_HPP
#define STAGETIMER_HPP

// System includes
//
#include <string>

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Timers/TimerMacro.h"

namespace GeoMHDiSCC {

   /**
    * @brief Implementation of a very simple stage timer
    */
   class StageTimer
   {
      public:
         /**
          * @brief New stage message
          *
          * @param msg Message to print
          */
         static void newStage(const std::string& msg);

         /**
          * @brief Stage completed message
          *
          * @param msg Message to print
          */
         static void completed(const std::string& msg);


         /**
          * @brief Stage message
          *
          * @param msg Message to print
          */
         static void msg(const std::string& msg);

         /**
          * @brief Start stage timing
          *
          * @param msg Message to print
          */
         static void start(const std::string& msg);

         /**
          * @brief End stage
          *
          * @param tabs Number of tab characters
          */
         static void done();

         /**
          * @brief Timer used for timings
          */
         static TimerMacro timer;
         
      protected:

      private:
         /**
          * @brief Constructor
          */
         StageTimer();

         /**
          * @brief Destructor
          */
         ~StageTimer();
   };

}

#endif // STAGETIMER_HPP
