/** \file IVectorEquation.hpp
 *  \brief Base for the implementation of a vector equation
 */

#ifndef IVECTOREQUATION_HPP
#define IVECTOREQUATION_HPP

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
#include "Enums/FieldComponents.hpp"
#include "Equations/EquationParameters.hpp"
#include "Equations/IEvolutionEquation.hpp"
#include "TypeSelectors/ScalarSelector.hpp"
#include "TypeSelectors/VariableSelector.hpp"

namespace GeoMHDiSCC {

   /**
    * \brief Base for the implementation of a vector equation
    */
   class IVectorEquation: public IEvolutionEquation
   {
      public:
         /**
          * \brief Simple constructor
          *
          * \param spEqParams Shared equation parameters
          */
         explicit IVectorEquation(SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~IVectorEquation();

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
          * @brief Compute the nonlinear interaction term
          *
          * @param rNLComp Nonlinear term component
          * @param name    ID of the physical vector component
          */
         virtual void computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) = 0;

         /**
          * @brief Compute the linear term
          *
          * @param rRHS RHS of timestepping equation
          * @param id   ID of the spectral vector component
          */
         virtual void computeLinear(Datatypes::SpectralScalarType& rRHS, FieldComponents::Spectral::Id id);

         /**
          * @brief Prepare the RHS for the timestep computation
          *
          * @param rhs    RHS of timestepping equation
          * @param compID  ID of the vector component
          */
         virtual void prepareTimestep(const Datatypes::SpectralScalarType& rhs, FieldComponents::Spectral::Id id) = 0;

      protected:
         /**
          * @brief Storage for the shared scalar variable
          */
         Datatypes::SharedVectorVariableType mspUnknown;

      private:
   };

   inline void IVectorEquation::setUnknown(Datatypes::SharedVectorVariableType spUnknown)
   {
      this->mspUnknown = spUnknown;
   }

   inline const Datatypes::VectorVariableType& IVectorEquation::unknown() const
   {
      return *this->mspUnknown;
   }

   inline Datatypes::VectorVariableType& IVectorEquation::rUnknown()
   {
      return *this->mspUnknown;
   }

   inline void IVectorEquation::computeLinear(Datatypes::SpectralScalarType& rRHS, FieldComponents::Spectral::Id id)
   {
   }

   /// Typedef for a shared IVectorEquation
   typedef SharedPtrMacro<IVectorEquation> SharedIVectorEquation;

}

#endif // IVECTOREQUATION_HPP
