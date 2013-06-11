/** \file TestTFTScalarTwo.hpp
 *  \brief Implementation of the second scalar equation to test TFT scheme
 */

#ifndef TESTTFTSCALARTWO_HPP
#define TESTTFTSCALARTWO_HPP

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
#include "Equations/Tests/ITestTFTScalarEquation.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   /**
    *  @brief Implementation of the second scalar equation to test TFT scheme
    */
   class TestTFTScalarTwo: public ITestTFTScalarEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams  Shared equation parameters
          */
         TestTFTScalarTwo(SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~TestTFTScalarTwo();

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

      private:
   };

}
}

#endif // TESTTFTSCALARTWO_HPP
