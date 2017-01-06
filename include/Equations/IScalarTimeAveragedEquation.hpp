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
#include "TypeSelectors/ScalarSelector.hpp"
#include "Equations/IScalarEquation.hpp"

namespace QuICC {

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

         /**
          * @brief Set the shared pointer to the unknown field
          *
          * This is required because the field are not initialised at creation time
          *
          * \param spUnknown Shared pointer to the unknown of the equation
          */
         virtual void setUnknown(Datatypes::SharedScalarVariableType spUnknown);
         
      protected:
         /**
          * @brief Update the stored value with the solver solution (read data)
          */
         virtual MHDFloat updateStoredSolution(const MHDFloat newData, FieldComponents::Spectral::Id compId, const int i, const int j, const int k);

         /**
          * @brief Update the stored value with the solver solution (complex data)
          */
         virtual MHDComplex updateStoredSolution(const MHDComplex newData, FieldComponents::Spectral::Id compId, const int i, const int j, const int k);

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
         SharedPtrMacro<Datatypes::SpectralScalarType> mTimeAvg;
   };

   /// Typedef for a shared IScalarTimeAveragedEquation
   typedef SharedPtrMacro<IScalarTimeAveragedEquation> SharedIScalarTimeAveragedEquation;
}
}

#endif // ISCALARTIMEAVERAGEDEQUATION_HPP
