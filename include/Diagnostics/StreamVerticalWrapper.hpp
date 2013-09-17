/**
 * @file StreamVerticalWrapper.hpp
 * @brief Streamfunction and vertical velocity wrapper into velocity field 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef STREAMVERTICALWRAPPER_HPP
#define STREAMVERTICALWRAPPER_HPP

// Configuration includes
//

// System includes
//
#include <string>

// External includes
//

// Project includes
//
#include "Diagnostics/IVelocityWrapper.hpp"
#include "TypeSelectors/VariableSelector.hpp"

namespace GeoMHDiSCC {

namespace Diagnostics {

   /**
    * @brief Streamfunction and vertical velocity wrapper into velocity field
    */
   class StreamVerticalWrapper: public IVelocityWrapper
   {
      public:
         /**
          * @brief Constructor
          */
         StreamVerticalWrapper(const Datatypes::SharedScalarVariableType spStream, const Datatypes::SharedScalarVariableType spVertical);

         /**
          * @brief Constructor
          */
         ~StreamVerticalWrapper();

         /**
          * @brief Get first velocity field component
          */
         virtual const Datatypes::PhysicalScalarType& one() const;

         /**
          * @brief Get second velocity field component
          */
         virtual const Datatypes::PhysicalScalarType& two() const;

         /**
          * @brief Get third velocity field component
          */
         virtual const Datatypes::PhysicalScalarType& three() const;

      protected:

      private:
         /**
          * @brief Shared streamfunction variable
          */
         Datatypes::SharedScalarVariableType mspStream;

         /**
          * @brief Shared vertical velocity variable
          */
         Datatypes::SharedScalarVariableType mspVertical;
   };

   /// Typedef for a shared StreamVerticalWrapper
   typedef SharedPtrMacro<StreamVerticalWrapper> SharedStreamVerticalWrapper;
}
}

#endif // STREAMVERTICALWRAPPER_HPP
