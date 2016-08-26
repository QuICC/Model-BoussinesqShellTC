/** 
 * @file IVectorTimeAveragedEquation.cpp
 * @brief Source of the base implementation of a vector time averaged equation
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "Equations/IVectorTimeAveragedEquation.hpp"

// Project includes
//
#include "Exceptions/Exception.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   IVectorTimeAveragedEquation::IVectorTimeAveragedEquation(SharedEquationParameters spEqParams)
      : IVectorEquation(spEqParams), mTimeFinished(true), mTimestep(-1.0)
   {
   }

   IVectorTimeAveragedEquation::~IVectorTimeAveragedEquation()
   {
   }

   void IVectorTimeAveragedEquation::setTime(const MHDFloat time, const bool finished)
   {
      this->mTimeFinished = finished;

      if(this->mTimeFinished)
      {
         this->mTimestep = time - this->time();
         EquationData::setTime(time, finished);
      }
   }

   MHDFloat IVectorTimeAveragedEquation::updateStoredSolution(const MHDFloat oldData, const MHDFloat newData)
   {
      // Only update mean on full timestep
      if(this->mTimeFinished)
      {
         return this->incrementTimeAverage(oldData, newData, this->time(), this->mTimestep);
      } else
      {
         return oldData;
      }
   }

   MHDComplex IVectorTimeAveragedEquation::updateStoredSolution(const MHDComplex oldData, const MHDComplex newData)
   {
      // Only update mean on full timestep
      if(this->mTimeFinished)
      {
         return this->incrementTimeAverage(oldData, newData, this->time(), this->mTimestep);
      } else
      {
         return oldData;
      }
   }
}
}
