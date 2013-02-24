/** \file VectorEquation.hpp
 *  \brief Base for the implementation of a vector equation
 */

#ifndef VECTOREQUATION_HPP
#define VECTOREQUATION_HPP

// Configuration includes
//
#include "SmartPointers/SharedPtrMacro.h"
#include "Simulation/PrepMacros/ScalarTypedefsMacro.h"
#include "Simulation/PrepMacros/VariableTypedefsMacro.h"

// System includes
//

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Enums/FieldComponents.h"
#include "Equations/EquationParameters.hpp"
#include "Equations/EvolutionEquation.hpp"

namespace GeoMHDiSCC {

   /**
    * \brief Base for the implementation of a vector equation
    */
   class VectorEquation: public EvolutionEquation
   {
      public:
         /**
          * \brief Simple constructor
          *
          * \param spEqParams Shared equation parameters
          */
         VectorEquation(SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~VectorEquation();

         /**
          * @brief Set the smart pointer to the unknown field
          *
          * This is required because the field are not initialised at creation time
          *
          * \param spUnknown Shared pointer to the unknown of the equation
          */
         void setUnknown(Code::SharedVectorVariable spUnknown);
         
         /**
          * @brief Get the unknown variable
          */
         const Code::VectorVariable& unknown() const;

         /**
          * @brief Set the unknown variable
          */
         Code::VectorVariable& rUnknown();

         /**
          * @brief Compute the nonlinear interaction term
          *
          * @param rNLComp Nonlinear term component
          * @param name    ID of the physical vector component
          */
         virtual void computeNonlinear(Code::PhysicalScalarType& rNLComp, FieldComponents::Physical::Component id) = 0;

         /**
          * @brief Compute the linear term
          *
          * @param rRHS RHS of timestepping equation
          * @param id   ID of the spectral vector component
          */
         virtual void computeLinear(Code::SpectralScalarType& rRHS, FieldComponents::Spectral::Component id);

         /**
          * @brief Prepare the RHS for the timestep computation
          *
          * @param rhs    RHS of timestepping equation
          * @param compID  ID of the vector component
          */
         virtual void prepareTimestep(const Code::SpectralScalarType& rhs, FieldComponents::Spectral::Component id) = 0;

      protected:
         /**
          * @brief Storage for the shared scalar variable
          */
         Code::SharedVectorVariable mspUnknown;

      private:
   };

   inline void VectorEquation::setUnknown(Code::SharedVectorVariable spUnknown)
   {
      this->mspUnknown = spUnknown;
   }

   inline const Code::VectorVariable& VectorEquation::unknown() const
   {
      return *this->mspUnknown;
   }

   inline Code::VectorVariable& VectorEquation::rUnknown()
   {
      return *this->mspUnknown;
   }

   inline void VectorEquation::computeLinear(Code::SpectralScalarType& rRHS, FieldComponents::Spectral::Component id)
   {
   }

   /// Typedef for a shared VectorEquation
   typedef SharedPtrMacro<VectorEquation> SharedVectorEquation;

}

#endif // VECTOREQUATION_HPP
