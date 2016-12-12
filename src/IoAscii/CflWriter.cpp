/** 
 * @file CflWriter.cpp 
 * @brief Source of the implementation of the CFL writer
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Configuration includes
//
#include "Framework/FrameworkMacro.h"

// System includes
//
#include <iomanip>

// External includes
//

// Class include
//
#include "IoAscii/CflWriter.hpp"

// Project includes
//
#include "IoAscii/CflTags.hpp"

namespace GeoMHDiSCC {

namespace IoAscii {

   CflWriter::CflWriter()
      : IAsciiEWriter(CflTags::NAME, CflTags::EXTENSION, CflTags::HEADER, CflTags::TYPE, CflTags::VERSION), mTime(-1.0), mTimestep(-1.0), mSteps(0.0), mChanged(false)
   {
   }

   CflWriter::~CflWriter()
   {
   }

   void CflWriter::setSimTime(const MHDFloat time, const MHDFloat timestep, const MHDFloat steps, const MHDFloat error)
   {
      this->mTime = time;

      this->mChanged = (timestep != this->mTimestep);
      this->mTimestep = timestep;

      this->mSteps = steps;

      this->mError = error;
   }

   void CflWriter::write()
   {
      // pre write
      this->preWrite();

      // Check if the workflow allows IO to be performed
      if(this->mChanged && FrameworkMacro::allowsIO())
      {
         this->mFile << std::setprecision(14) << this->mTime << "\t" << this->mTimestep << "\t" << this->mSteps << "\t" << this->mError << std::endl;
      }

      // post write
      this->postWrite();
   }

}
}
