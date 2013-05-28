/** \file TestTFTScalarThree.hpp
 *  \brief Implementation of the third scalar equation to test TFT scheme
 */

#ifndef TESTTFTSCALARTHREE_HPP
#define TESTTFTSCALARTHREE_HPP

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
    *  @brief Implementation of the third scalar equation to test TFT scheme
    */
   class TestTFTScalarThree: public ITestTFTScalarEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams  Shared equation parameters
          */
         TestTFTScalarThree(SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~TestTFTScalarThree();

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

#endif // TESTTFTSCALARTHREE_HPP