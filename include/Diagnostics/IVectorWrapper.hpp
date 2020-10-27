/**
 * @file IVectorWrapper.hpp
 * @brief Interface to the vector wrappers used by the diagnostics 
 */

#ifndef IVECTORWRAPPER_HPP
#define IVECTORWRAPPER_HPP

// Configuration includes
//
#include "SmartPointers/SharedPtrMacro.h"

// System includes
//
#include <string>

// External includes
//

// Project includes
//
#include "TypeSelectors/VariableSelector.hpp"

namespace QuICC {

namespace Diagnostics {

   /**
    * @brief Interface to the vector field wrappers used by the diagnostics
    */
   class IVectorWrapper
   {
      public:
         /**
          * @brief Constructor
          */
         IVectorWrapper();

         /**
          * @brief Constructor
          */
         ~IVectorWrapper();

         /**
          * @brief Get first vector field component
          */
         virtual const Datatypes::PhysicalScalarType& one() const = 0;

         /**
          * @brief Get second vector field component
          */
         virtual const Datatypes::PhysicalScalarType& two() const = 0;

         /**
          * @brief Get third vector field component
          */
         virtual const Datatypes::PhysicalScalarType& three() const = 0;

         /**
          * @brief Get Resolution
          */
         virtual const SharedResolution spRes() const = 0;

      protected:

      private:
   };

   /// Typedef for a shared IVectorWrapper
   typedef SharedPtrMacro<IVectorWrapper> SharedIVectorWrapper;
}
}

#endif // IVECTORWRAPPER_HPP
