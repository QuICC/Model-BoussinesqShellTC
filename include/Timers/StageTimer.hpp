/**
 * @file StageTimer.hpp
 * @brief Implementation of a very simple stage timer 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef STAGETIMER_HPP
#define STAGETIMER_HPP

// Configuration includes
//
#include "Timers/TimerMacro.h"

// System includes
//
#include <string>

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"

namespace GeoMHDiSCC {

   /**
    * @brief Implementation of a very simple stage timer
    */
   class StageTimer
   {
      public:
         /**
          * @brief Constructor
          */
         explicit StageTimer();

         /**
          * @brief Destructor
          */
         ~StageTimer();

         /**
          * @brief New stage message
          *
          * @param msg Message to print
          */
         static void stage(const std::string& msg);

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
         static void msg(const std::string& msg, const int space = 8);

         /**
          * @brief Start stage timing
          *
          * @param msg Message to print
          */
         void start(const std::string& msg, const int level = 0);

         /**
          * @brief End stage
          *
          * @param tabs Number of tab characters
          */
         void done();
         
      protected:

      private:
         /**
          * @brief Level
          */
         int mLevel;

         /**
          * @brief The actual timer
          */
         TimerMacro mTimer;
   };

}

#endif // STAGETIMER_HPP
