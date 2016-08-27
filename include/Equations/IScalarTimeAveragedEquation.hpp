/**
 * @file IScalarTimeAveragedEquation.hpp
 * @brief Base for the implementation of a time averaged scalar equation  
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef ISCALARTIMEAVERAGEDEQUATION_HPP
#define ISCALARTIMEAVERAGEDEQUATION_HPP

// Configuration includes
//
#include "SmartPointers/SharedPtrMacro.h"

// System includes
//

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Equations/IScalarEquation.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   /**
    * @brief Base for the implementation of a time averaged scalar equation
    */
   class IScalarTimeAveragedEquation: public IScalarEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * \param spEqParams Shared equation parameters
          */
         explicit IScalarTimeAveragedEquation(SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~IScalarTimeAveragedEquation();

         /**
          * @brief Current simulation time update
          */
         virtual void setTime(const MHDFloat time, const bool finished);
         
      protected:
         /**
          * @brief Update the stored value with the solver solution (read data)
          */
         virtual MHDFloat updateStoredSolution(const MHDFloat oldData, const MHDFloat newData);

         /**
          * @brief Update the stored value with the solver solution (complex data)
          */
         virtual MHDComplex updateStoredSolution(const MHDComplex oldData, const MHDComplex newData);

      private:
         /**
          * @brief Multi staged timestep finished?
          */
         bool mTimeFinished;

         /**
          * @brief Timestep since last update
          */
         MHDFloat mTimestep;

         /**
          * @brief Storage of the previous values
          */
         SharedSpectralScalarType mTimeAvg;

   };

   /// Typedef for a shared IScalarTimeAveragedEquation
   typedef SharedPtrMacro<IScalarTimeAveragedEquation> SharedIScalarTimeAveragedEquation;
}
}

#endif // ISCALARTIMEAVERAGEDEQUATION_HPP
