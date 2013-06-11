/** \file IScalarDEquation.hpp
 *  \brief Base for the implementation of a scalar diagnostic equation
 */

#ifndef ISCALARDEQUATION_HPP
#define ISCALARDEQUATION_HPP

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
#include "Equations/EquationParameters.hpp"
#include "Equations/IDiagnosticEquation.hpp"
#include "TypeSelectors/ScalarSelector.hpp"
#include "TypeSelectors/VariableSelector.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   /**
    * \brief Base for the implementation of a scalar diagnostic equation
    */
   class IScalarDEquation: public IDiagnosticEquation
   {
      public:
         /**
          * \brief Simple constructor
          *
          * \param spEqParams Shared equation parameters
          */
         explicit IScalarDEquation(SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~IScalarDEquation();

         /**
          * @brief Set the shared pointer to the unknown field
          *
          * This is required because the field are not initialised at creation time
          *
          * \param spUnknown Shared pointer to the unknown of the equation
          */
         virtual void setUnknown(Datatypes::SharedScalarVariableType spUnknown);

         /**
          * @brief Get the unknown variable
          */
         const Datatypes::ScalarVariableType& unknown() const;

         /**
          * @brief Set the unknown variable
          */
         Datatypes::ScalarVariableType& rUnknown();
         
      protected:

      private:
         /**
          * @brief Storage for the shared scalar variable
          */
         Datatypes::SharedScalarVariableType mspUnknown;
   };

   /// Typedef for a shared IScalarDEquation
   typedef SharedPtrMacro<IScalarDEquation> SharedIScalarDEquation;

}
}

#endif // ISCALARDEQUATION_HPP
