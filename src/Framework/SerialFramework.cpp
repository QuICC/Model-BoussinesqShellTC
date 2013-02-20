/** \file SerialFramework.cpp
 *  \brief Source of implementation of serial framework
 */

// Configuration includes
//
#include "Profiler/ProfilerMacro.h"

// System includes
//

// External includes
//

// Class include
//
#include "Base/Framework/SerialFramework.hpp"

// Project includes
//
#include "Base/Precision.hpp"

namespace GeoMHDiSCC {

   void SerialFramework::init()
   {
      // Initialise the profiler if needed
      ProfilerMacro_init();

      // Initialise the precision framework
      Precision::init();
   }

   void SerialFramework::setup(const int nCpu)
   {
      // Set the number of CPUs
      SerialFramework::mNCpu = nCpu;

      // Set ID
      SerialFramework::mCpuId = 0;

      // Check that the framework setup is right
      SerialFramework::checkFramework(1);
   }

   void SerialFramework::synchronize()
   {
      // Nothing to be done in serial framework
   }

   void SerialFramework::finalize()
   {
      // Nothing to be done in serial framework
   }

   SerialFramework::SerialFramework()
   {
   }

}
