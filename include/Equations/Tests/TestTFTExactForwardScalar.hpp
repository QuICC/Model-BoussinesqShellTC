/**
 * @file TestTFTExactForwardScalar.hpp
 * @brief Implementation of a test equation for the TFT scheme with exact known physical space scalar solution
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef TESTTFTEXACTFORWARDSCALAR_HPP
#define TESTTFTEXACTFORWARDSCALAR_HPP

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

namespace GeoMHDiSCC {

namespace Equations {

   /**
    * @brief Implementation of a test equation for the TFT scheme with exact known physical space scalar solution
    */
   class TestTFTExactForwardScalar: public IScalarEquation
   {
      public:
         /**
          * @brief Constructor
          */
         TestTFTExactForwardScalar();

         /**
          * @brief Destructor
          */
         virtual ~TestTFTExactForwardScalar();

         /**
          * @brief Set the unknown name and requirements
          */
         void setIdentity(const PhysicalNames::Id name);

         /**
          * @brief Set the exact solution as nonlinear interaction term
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

   /// Typedef for a shared TestTFTExactForwardScalar
   typedef SharedPtrMacro<TestTFTExactForwardScalar> SharedTestTFTExactForwardScalar;

}
}

#endif // TESTTFTEXACTFORWARDSCALAR_HPP
