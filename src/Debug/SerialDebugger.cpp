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

   void SerialDebugger::enter(const std::string& msg)
   {
      std::cerr << "ENTERING: " << msg << std::endl;
   }

   void SerialDebugger::leave(const std::string& msg)
   {
      std::cerr << "LEAVING: " << msg << std::endl;
   }

}
}
