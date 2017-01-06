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

namespace QuICC {

namespace Equations {

   IVectorTimeAveragedEquation::IVectorTimeAveragedEquation(SharedEquationParameters spEqParams)
      : IVectorEquation(spEqParams), mTimeFinished(true), mTimestep(-1.0)
   {
      EquationData::setTime(-4242.0, false);
   }

   IVectorTimeAveragedEquation::~IVectorTimeAveragedEquation()
   {
   }

   void IVectorTimeAveragedEquation::setTime(const MHDFloat time, const bool finished)
   {
      if(this->time() == -4242.0)
      {
         for(SpectralComponent_iterator it = this->spectralRange().first; it != this->spectralRange().second; it++)
         {
            this->mTimeAvg->rComp(*it).setData(this->unknown().dom(0).perturbation().comp(*it).data());
         }
      }

      this->mTimeFinished = finished;

      if(this->mTimeFinished)
      {
         this->mTimestep = time - this->time();
         EquationData::setTime(time, finished);
      }
   }

   void IVectorTimeAveragedEquation::setUnknown(Datatypes::SharedVectorVariableType spUnknown)
   {
      IVectorEquation::setUnknown(spUnknown);

      this->mTimeAvg = SharedPtrMacro<Datatypes::VectorField<Datatypes::SpectralScalarType,FieldComponents::Spectral::Id> >(new Datatypes::VectorField<Datatypes::SpectralScalarType,FieldComponents::Spectral::Id>(this->unknown().dom(0).perturbation()));
   }

   MHDFloat IVectorTimeAveragedEquation::updateStoredSolution(const MHDFloat newData, FieldComponents::Spectral::Id compId, const int i, const int j, const int k)
   {
      // Only update mean on full timestep
      if(this->mTimeFinished)
      {
         MHDFloat val = incrementTimeAverage(this->mTimeAvg->comp(compId).point(i, j, k), newData, this->time(), this->mTimestep);
         this->mTimeAvg->rComp(compId).setPoint(val, i, j, k);
         return val;
      } else
      {
         return noupdateTimeAverage(this->mTimeAvg->comp(compId).point(i, j, k), newData);
      }
   }

   MHDComplex IVectorTimeAveragedEquation::updateStoredSolution(const MHDComplex newData, FieldComponents::Spectral::Id compId, const int i, const int j, const int k)
   {
      // Only update mean on full timestep
      if(this->mTimeFinished)
      {
         MHDComplex val = incrementTimeAverage(this->mTimeAvg->comp(compId).point(i, j, k), newData, this->time(), this->mTimestep);
         this->mTimeAvg->rComp(compId).setPoint(val, i, j, k);
         return val;
      } else
      {
         return noupdateTimeAverage(this->mTimeAvg->comp(compId).point(i, j, k), newData);
      }
   }
}
}
