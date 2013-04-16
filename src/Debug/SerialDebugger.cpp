/** \file SerialDebugger.cpp
 *  \brief Source of the serial debugger implementation
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

   void SerialDebugger::enter(const std::string& msg, const int tabs)
   {
      for(int i = 0; i < tabs; i++)
      {
         std::cerr << "\t";
      }
      std::cerr << "ENTERING: " << msg << std::endl;
   }

   void SerialDebugger::leave(const std::string& msg, const int tabs)
   {
      for(int i = 0; i < tabs; i++)
      {
         std::cerr << "\t";
      }
      std::cerr << "LEAVING: " << msg << std::endl;
   }

}
}
