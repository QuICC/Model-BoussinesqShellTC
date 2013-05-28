/** \file TestTFTScalarOne.hpp
 *  \brief Implementation of the first scalar equation to test TFT scheme
 */

#ifndef TESTTFTSCALARONE_HPP
#define TESTTFTSCALARONE_HPP

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
    *  @brief Implementation of the first scalar equation to test TFT scheme
    */
   class TestTFTScalarOne: public ITestTFTScalarEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams  Shared equation parameters
          */
         TestTFTScalarOne(SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~TestTFTScalarOne();

         /**
          * @brief Compute the nonlinear interaction term
          *
          * @param rNLComp Nonlinear term component
          */
         virtual void computeNonlinear(Datatypes::PhysicalScalarType& rNLComp) const;

      protected:
         /**
          * @brief Set variable requirements
          */
         virtual void setRequirements();

      private:
   };

}
}

#endif // TESTTFTSCALARONE_HPP
