/** \file FrameworkBase.cpp
 *  \brief Source of base of implementation of framework
 */

// System includes
//

// External includes
//

// Class include
//
#include "Framework/FrameworkBase.hpp"

// Project includes
//
#include "Exceptions/Exception.hpp"

namespace GeoMHDiSCC {

   const int FrameworkBase::IO_RANK = 0;

   int FrameworkBase::mNCpu = -1;

   int FrameworkBase::mCpuId = -1;

   void FrameworkBase::checkFramework(const int nCpu)
   {
      // Check that the number of CPUs was set
      if(FrameworkBase::nCpu() < 1)
      {
         throw Exception("Framework contains a negative number of CPUs!");
      }

      // Check that the CPU ID was set
      if(FrameworkBase::id() < 0)
      {
         throw Exception("Framework CPU ID was not set!");
      }

      // Check compatibility between requested cores and setup cores
      if(nCpu != FrameworkBase::nCpu())
      {
         throw Exception("Framework and parameters have conflicting number of CPUs");
      }
   }

   bool FrameworkBase::allowsIO()
   {
      // Only IO_RANK is allowed to do (direct) IO
      return (FrameworkBase::id() == FrameworkBase::IO_RANK);
   }

   FrameworkBase::FrameworkBase()
   {
   }

   FrameworkBase::~FrameworkBase()
   {
   }
}
