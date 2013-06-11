/** \file IVectorDEquation.hpp
 *  \brief Base for the implementation of a vector diagnostic equation
 */

#ifndef IVECTORDEQUATION_HPP
#define IVECTORDEQUATION_HPP

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
#include "Enums/FieldIds.hpp"
#include "Equations/EquationParameters.hpp"
#include "Equations/IDiagnosticEquation.hpp"
#include "TypeSelectors/ScalarSelector.hpp"
#include "TypeSelectors/VariableSelector.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   /**
    * \brief Base for the implementation of a vector diagnostic equation
    */
   class IVectorDEquation: public IDiagnosticEquation
   {
      public:
         /**
          * \brief Simple constructor
          *
          * \param spEqParams Shared equation parameters
          */
         explicit IVectorDEquation(SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~IVectorDEquation();

         /**
          * @brief Set the smart pointer to the unknown field
          *
          * This is required because the field are not initialised at creation time
          *
          * \param spUnknown Shared pointer to the unknown of the equation
          */
         void setUnknown(Datatypes::SharedVectorVariableType spUnknown);
         
         /**
          * @brief Get the unknown variable
          */
         const Datatypes::VectorVariableType& unknown() const;

         /**
          * @brief Set the unknown variable
          */
         Datatypes::VectorVariableType& rUnknown();

      protected:

      private:
         /**
          * @brief Storage for the shared scalar variable
          */
         Datatypes::SharedVectorVariableType mspUnknown;
   };

   /// Typedef for a shared IVectorDEquation
   typedef SharedPtrMacro<IVectorDEquation> SharedIVectorDEquation;

}
}

#endif // IVECTORDEQUATION_HPP
