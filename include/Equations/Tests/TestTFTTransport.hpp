/** \file TestTFTTransport.hpp
 *  \brief Implementation of the transport equation to test TFT scheme
 */

#ifndef TESTTFTTRANSPORT_HPP
#define TESTTFTTRANSPORT_HPP

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
    *  @brief Implementation of the transport equation to test TFT scheme
    */
   class TestTFTTransport: public ITestTFTScalarEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams  Shared equation parameters
          */
         TestTFTTransport(SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~TestTFTTransport();

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

#endif // TESTTFTTRANSPORT_HPP
