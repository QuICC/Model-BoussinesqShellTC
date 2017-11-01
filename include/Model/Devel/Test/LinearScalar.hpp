/**
 * @file LinearScalar.hpp
 * @brief Implementation of the scalar linear test equation
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef TESTLINEARSCALAR_HPP
#define TESTLINEARSCALAR_HPP

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
#include "Equations/IScalarEquation.hpp"

namespace QuICC {

namespace Equations {

namespace Test {

   /**
    * @brief Implementation of the scalar linear test equation
    */
   class LinearScalar: public IScalarEquation
   {
      public:
         /**
          * @brief Constructor
          *
          * @param spEqParams  Shared equation parameters
          */
         LinearScalar(SharedEquationParameters spEqParams);

         /**
          * @brief Destructor
          */
         virtual ~LinearScalar();

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

   /// Typedef for a shared LinearScalar
   typedef SharedPtrMacro<LinearScalar> SharedLinearScalar;

}
}
}

#endif // TESTLINEARSCALAR_HPP