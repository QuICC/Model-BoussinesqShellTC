/** \file IVectorPEquation.hpp
 *  \brief Base for the implementation of a vector prognostic equation
 */

#ifndef IVECTORPEQUATION_HPP
#define IVECTORPEQUATION_HPP

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
#include "Equations/IPrognosticEquation.hpp"
#include "TypeSelectors/ScalarSelector.hpp"
#include "TypeSelectors/VariableSelector.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   /**
    * \brief Base for the implementation of a vector prognostic equation
    */
   class IVectorPEquation: public IPrognosticEquation
   {
      public:
         /**
          * \brief Simple constructor
          *
          * \param spEqParams Shared equation parameters
          */
         explicit IVectorPEquation(SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~IVectorPEquation();

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

         /**
          * @brief Prepare the RHS for the timestep computation
          *
          * @param rhs  RHS of timestepping equation
          * @param id   ID of the vector component
          */
         virtual void prepareTimestep(const Datatypes::SpectralScalarType& rhs, FieldComponents::Spectral::Id id);

      protected:

      private:
         /**
          * @brief Storage for the shared scalar variable
          */
         Datatypes::SharedVectorVariableType mspUnknown;
   };

   /// Typedef for a shared IVectorPEquation
   typedef SharedPtrMacro<IVectorPEquation> SharedIVectorPEquation;

}
}

#endif // IVECTORPEQUATION_HPP
