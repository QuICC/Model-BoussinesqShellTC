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

namespace GeoMHDiSCC {

namespace Equations {

   IScalarTimeAveragedEquation::IScalarTimeAveragedEquation(SharedEquationParameters spEqParams)
      : IScalarEquation(spEqParams), mTimeFinished(false), mTimestep(-1.0)
   {
      EquationData::setTime(-4242.0, false);
   }

   IScalarTimeAveragedEquation::~IScalarTimeAveragedEquation()
   {
   }

   void IScalarTimeAveragedEquation::setTime(const MHDFloat time, const bool finished)
   {
      if(this->time() == -4242.0)
      {
         this->mTimeAvg->setData(this->unknown().dom(0).perturbation().data());
         EquationData::setTime(time, finished);
      }

      this->mTimeFinished = finished;

      if(this->mTimeFinished)
      {
         this->mTimestep = time - this->time();
         EquationData::setTime(time, finished);
      }
   }

   void IScalarTimeAveragedEquation::setUnknown(Datatypes::SharedScalarVariableType spUnknown)
   {
      IScalarEquation::setUnknown(spUnknown);

      this->mTimeAvg = SharedPtrMacro<Datatypes::SpectralScalarType>(new Datatypes::SpectralScalarType(this->unknown().dom(0).perturbation()));
   }

   MHDFloat IScalarTimeAveragedEquation::updateStoredSolution(const MHDFloat newData, FieldComponents::Spectral::Id compId, const int i, const int j, const int k)
   {
      // Only update mean on full timestep
      if(this->mTimeFinished)
      {
         MHDFloat val = incrementTimeAverage(this->mTimeAvg->point(i, j, k), newData, this->time(), this->mTimestep);
         this->mTimeAvg->setPoint(val, i, j, k);
         return val; 
      } else
      {
         return noupdateTimeAverage(this->mTimeAvg->point(i,j,k), newData);
      }
   }

   MHDComplex IScalarTimeAveragedEquation::updateStoredSolution(const MHDComplex newData, FieldComponents::Spectral::Id compId, const int i, const int j, const int k)
   {
      // Only update mean on full timestep
      if(this->mTimeFinished)
      {
         MHDComplex val = incrementTimeAverage(this->mTimeAvg->point(i, j, k), newData, this->time(), this->mTimestep);
         this->mTimeAvg->setPoint(val, i, j, k);
         return val; 
      } else
      {
         return noupdateTimeAverage(this->mTimeAvg->point(i,j,k), newData);
      }
   }

}
}
