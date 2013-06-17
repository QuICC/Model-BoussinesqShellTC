/** \file FieldVisualizer.hpp
 *  \brief Implementation of the basic field visualizer
 */

#ifndef FIELDVISUALIZER_HPP
#define FIELDVISUALIZER_HPP

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
    * @brief Implementation of thecw basic field visualizer
    */
   class FieldVisualizer: public IScalarEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams Shared equation parameters
          * @paarm name       Name of the field
          */
         FieldVisualizer(SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~FieldVisualizer();

         /**
          * @brief Set the unknown name and requirements
          */
         void setIdentity(const PhysicalNames::Id name);

         /**
          * @brief Set which fields to output
          */
         void setFields(const bool viewField, const bool viewGradient);

      protected:
         /**
          * @brief Set variable requirements
          */
         virtual void setRequirements();

         /**
          * @brief Set coupling information
          */
         virtual void setCoupling();

      private:
         /**
          * @brief Storage for output field flag
          */
         bool mViewField;

         /**
          * @brief Storage for output gradient flag
          */
         bool mViewGradient;
   };

   /// Typedef for a shared FieldVisualizer
   typedef SharedPtrMacro<FieldVisualizer> SharedFieldVisualizer;

}
}

#endif // FIELDVISUALIZER_HPP
