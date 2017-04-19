/**
 * @file TimeAvgScalar.hpp
 * @brief Implementation of the time averaged equation for the timestepper test in TFF scheme
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef QUICC_EQUATIONS_TEST_TFF_TIMESTEPPER_TIMEAVGSCALAR_HPP
#define QUICC_EQUATIONS_TEST_TFF_TIMESTEPPER_TIMEAVGSCALAR_HPP

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

namespace Test {

namespace TFF {

namespace Timestepper {

   /**
    * @brief Implementation of the time averaged test equation for the timestepper test in TFF scheme
    */
   class TimeAvgScalar: public IScalarTimeAveragedEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams    Shared equation parameters
          */
         TimeAvgScalar(SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~TimeAvgScalar();

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
}
}
}

#endif // QUICC_EQUATIONS_TEST_TFF_TIMESTEPPER_TIMEAVGSCALAR_HPP
