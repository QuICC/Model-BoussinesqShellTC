/**
 * @file IVectorTimeAveragedEquation.hpp
 * @brief Base for the implementation of a vector time averaged equation 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef IVECTORTIMEAVERAGEDEQUATION_HPP
#define IVECTORTIMEAVERAGEDEQUATION_HPP

// Configuration includes
//
#include "SmartPointers/SharedPtrMacro.h"

// System includes
//
#include <vector>

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Equations/EquationParameters.hpp"
#include "Equations/IVectorEquation.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   /**
    * @brief Base for the implementation of a vector time averaged equation
    */
   class IVectorTimeAveragedEquation: public IVectorEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * \param spEqParams Shared equation parameters
          */
         explicit IVectorTimeAveragedEquation(SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~IVectorTimeAveragedEquation();

         /**
          * @brief Current simulation time update
          */
         virtual void setTime(const MHDFloat time, const bool finished);

      protected:
         /**
          * @brief Update the stored value with the solver solution (real data)
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
   };

   /// Typedef for a shared IVectorTimeAveragedEquation
   typedef SharedPtrMacro<IVectorTimeAveragedEquation> SharedIVectorTimeAveragedEquation;
}
}

#endif // IVECTORTIMEAVERAGEDEQUATION_HPP
