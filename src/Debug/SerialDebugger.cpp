/** 
 * @file SerialDebugger.cpp
 * @brief Source of the serial debugger implementation
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 * @version 0.9.0
 * @date 2013-09-11
 */

// System includes
//
#include <iostream>

// External includes
//

// Class include
//
#include "Debug/SerialDebugger.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Debug {

   TimerMacro SerialDebugger::timer = TimerMacro();

   void SerialDebugger::msg(const std::string& msg, const int tabs)
   {
      for(int i = 0; i < tabs; i++)
      {
         std::cerr << "\t";
      }
      std::cerr << "DEBUG " << msg << std::endl;
   }

   void SerialDebugger::enter(const std::string& msg, const int tabs)
   {
      SerialDebugger::msg("ENTERING: " + msg, tabs);
   }

   void SerialDebugger::leave(const std::string& msg, const int tabs)
   {
      SerialDebugger::msg("LEAVING: " + msg, tabs);
   }

   void SerialDebugger::showValue(const std::string& msg, const int tabs, const MHDFloat value)
   {
      for(int i = 0; i < tabs; i++)
      {
         std::cerr << "\t";
      }
      std::cerr << "DEBUG " << msg << value << std::endl;
   }

}
}
