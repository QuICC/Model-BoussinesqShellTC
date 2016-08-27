/** 
 * @file IScalarTimeAveragedEquation.cpp
 * @brief Source of the base implementation of a scalar time averaged equation
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
#include "Equations/IScalarTimeAveragedEquation.hpp"

// Project includes
//

#include <iostream>
namespace GeoMHDiSCC {

namespace Equations {

   IScalarTimeAveragedEquation::IScalarTimeAveragedEquation(SharedEquationParameters spEqParams)
      : IScalarEquation(spEqParams), mTimeFinished(true), mTimestep(-1.0)
   {
   }

   IScalarTimeAveragedEquation::~IScalarTimeAveragedEquation()
   {
   }

   void IScalarTimeAveragedEquation::setTime(const MHDFloat time, const bool finished)
   {
      this->mTimeFinished = finished;

      if(this->mTimeFinished)
      {
         this->mTimestep = time - this->time();
         EquationData::setTime(time, finished);
      }
   }

   MHDFloat IScalarTimeAveragedEquation::updateStoredSolution(const MHDFloat oldData, const MHDFloat newData)
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

   MHDComplex IScalarTimeAveragedEquation::updateStoredSolution(const MHDComplex oldData, const MHDComplex newData)
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
