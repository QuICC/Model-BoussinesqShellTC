/**
 * @file TestLinearVector.hpp
 * @brief Implementation of the vector linear test equation
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef TESTLINEARVECTOR_HPP
#define TESTLINEARVECTOR_HPP

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
#include "Equations/IVectorEquation.hpp"

namespace QuICC {

namespace Equations {

   /**
    * @brief Implementation of the vector linear test equation
    */
   class TestLinearVector: public IVectorEquation
   {
      public:
         /**
          * @brief Constructor
          *
          * @param spEqParams  Shared equation parameters
          */
         TestLinearVector(SharedEquationParameters spEqParams);

         /**
          * @brief Destructor
          */
         virtual ~TestLinearVector();

         /**
          * @brief Set the unknown name and requirements
          */
         void setIdentity(const PhysicalNames::Id name);

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

   /// Typedef for a shared TestLinearVector
   typedef SharedPtrMacro<TestLinearVector> SharedTestLinearVector;

}
}

#endif // TESTLINEARVECTOR_HPP
