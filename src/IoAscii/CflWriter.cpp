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
#include <sstream>

// External includes
//

// Class include
//
#include "IoAscii/CflWriter.hpp"

// Project includes
//
#include "IoAscii/CflTags.hpp"

namespace QuICC {

namespace IoAscii {

   CflWriter::CflWriter()
      : IAsciiWriter(CflTags::NAME, CflTags::EXTENSION, CflTags::HEADER, CflTags::TYPE, CflTags::VERSION, IAsciiWriter::EXTEND), mTime(-1.0), mTimestep(-1.0), mSteps(0.0), mChanged(false), mDt(2,1), mcIoHigh(14), mcIoLow(7), mcIoExpW(5), mcIoFixW(2), mcIoIntW(7), mNeedFancy(true)
   {
      mDt.setConstant(-1.0);
   }

   CflWriter::~CflWriter()
   {
   }

   void CflWriter::setSimTime(const MHDFloat time, const Matrix& dt, const MHDFloat steps)
   {
      this->mTime = time;

      this->mChanged = (dt(0,0) != this->mDt(0,0) || static_cast<long int>(this->mSteps) % 100 == 0);
      this->mDt = dt;

      this->mSteps = steps;
   }

   void CflWriter::write()
   {
      this->fancyHeader();

      // pre write
      this->preWrite();

      int ioW = this->mcIoLow+this->mcIoExpW;
      int ioWH = this->mcIoHigh+this->mcIoExpW;
      int ioWf = this->mcIoLow+this->mcIoFixW;
      int ioWHf = this->mcIoHigh+this->mcIoFixW;

      // Check if the workflow allows IO to be performed
      if(this->mChanged && FrameworkMacro::allowsIO())
      {
         this->mFile << std::setprecision(this->mcIoHigh);
         this->mFile << std::left << std::setw(ioWH) << std::scientific << this->mTime << "\t";

         this->mFile << std::left << std::setw(ioWH) << std::scientific << this->mDt(0,0) << "\t";
         for(int j = 1; j < this->mDt.rows(); j++)
         {
            if(this->mDt(j,0) == -1)
            {
               this->mFile << std::left << std::setfill(' ') << std::setw(ioWHf) << "global" << "\t";
            } else if(this->mDt(j,0) == -100)
            {
               this->mFile << std::left << std::setfill(' ') << std::setw(ioWHf) << "fixed" << "\t";
            } else if(this->mDt(j,0) == -101)
            {
               this->mFile << std::left << std::setfill(' ') << std::setw(ioWHf) << "max" << "\t";
            } else if(this->mDt(j,0) == -102)
            {
               this->mFile << std::left << std::setfill(' ') << std::setw(ioWHf) << "min" << "\t";
            } else
            {
               this->mFile << std::left << std::setw(ioWHf) << std::fixed << this->mDt(j,0) << "\t";
            }
         }
         this->mFile << std::setprecision(this->mcIoLow);
         for(int i = 1; i < this->mDt.cols(); i++)
         {
            this->mFile << std::left << std::setw(ioW) << std::scientific << this->mDt(0,i) << "\t";
            for(int j = 1; j < this->mDt.rows(); j++)
            {
               if(this->mDt(j,i) == -1)
               {
                  this->mFile << std::left << std::setfill(' ') << std::setw(ioW) << "global" << "\t";
               } else
               {
                  this->mFile << std::left << std::setw(ioWf) << std::fixed << this->mDt(j,i) << "\t";
               }
            }
         }
         this->mFile << std::defaultfloat;
         this->mFile << std::left << std::setw(this->mcIoIntW) << this->mSteps << std::endl;
      }

      // post write
      this->postWrite();
   }

   void CflWriter::fancyHeader()
   {
      if(this->mNeedFancy)
      {
         int ioW = this->mcIoLow+this->mcIoExpW;
         int ioWH = this->mcIoHigh+this->mcIoExpW;
         int ioWf = this->mcIoLow+this->mcIoFixW;
         int ioWHf = this->mcIoHigh+this->mcIoFixW;

         this->mFile << std::left << std::setfill(' ') << std::setw(ioWH) << "# time" << "\t";
         this->mFile << std::setw(ioWH) << "dt" << "\t";
         this->mFile << std::setw(ioWHf) << "loc" << "\t";

         std::stringstream ss;
         for(int i = 1; i < this->mDt.cols(); i++)
         {
            std::stringstream ss;
            ss << "dt_" << i;
            this->mFile << std::setw(ioW) << ss.str() << "\t";
            ss.str("");
            ss << "loc_" << i;
            this->mFile << std::setw(ioWf) << ss.str() << "\t";
         }

         this->mFile << std::setw(this->mcIoIntW) << "steps";
         this->mFile << std::endl;

         this->mNeedFancy = false;
      }
   }

}
}
