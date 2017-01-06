/**
 * @file IVelocityWrapper.hpp
 * @brief Interface to the velocity wrappers used by the diagnostics 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef IVELOCITYWRAPPER_HPP
#define IVELOCITYWRAPPER_HPP

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
    * @brief Interface to the velocity wrappers used by the diagnostics
    */
   class IVelocityWrapper
   {
      public:
         /**
          * @brief Constructor
          */
         IVelocityWrapper();

         /**
          * @brief Constructor
          */
         ~IVelocityWrapper();

         /**
          * @brief Get first velocity field component
          */
         virtual const Datatypes::PhysicalScalarType& one() const = 0;

         /**
          * @brief Get second velocity field component
          */
         virtual const Datatypes::PhysicalScalarType& two() const = 0;

         /**
          * @brief Get third velocity field component
          */
         virtual const Datatypes::PhysicalScalarType& three() const = 0;

         /**
          * @brief Get Resolution
          */
         virtual const SharedResolution spRes() const = 0;

      protected:

      private:
   };

   /// Typedef for a shared IVelocityWrapper
   typedef SharedPtrMacro<IVelocityWrapper> SharedIVelocityWrapper;
}
}

#endif // IVELOCITYWRAPPER_HPP
