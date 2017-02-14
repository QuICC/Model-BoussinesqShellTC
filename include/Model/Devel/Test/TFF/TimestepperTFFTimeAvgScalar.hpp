/**
 * @file TimestepperTFFTimeAvgScalar.hpp
 * @brief Implementation of the time averaged equation for the timestepper test in TFF scheme
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef TIMESTEPPERTFFTIMEAVGSCALAR_HPP
#define TIMESTEPPERTFFTIMEAVGSCALAR_HPP

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
#include "Equations/IScalarTimeAveragedEquation.hpp"

namespace QuICC {

namespace Equations {

   /**
    * @brief Implementation of the time averaged test equation for the timestepper test in TFF scheme
    */
   class TimestepperTFFTimeAvgScalar: public IScalarTimeAveragedEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams    Shared equation parameters
          */
         TimestepperTFFTimeAvgScalar(SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~TimestepperTFFTimeAvgScalar();

         /**
          * @brief Compute the nonlinear interaction term
          *
          * @param rNLComp Nonlinear term component
          * @param id      ID of the component (allows for a more general implementation)
          */
         virtual void computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const;

      protected:
         /**
          * @brief Set variable requirements
          */
         virtual void setRequirements();

         /**
          * @brief Set the equation coupling information
          */
         virtual void setCoupling();

      private:
   };

}
}

#endif // TIMESTEPPERTFFTIMEAVGSCALAR_HPP
