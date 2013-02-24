/** \file ScalarEquation.hpp
 *  \brief Base for the implementation of a scalar equation
 */

#ifndef SCALAREQUATION_HPP
#define SCALAREQUATION_HPP

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
#include "Equations/EquationParameters.hpp"
#include "Equations/EvolutionEquation.hpp"

namespace GeoMHDiSCC {

   /**
    * \brief Base for the implementation of a scalar equation
    */
   class ScalarEquation: public EvolutionEquation
   {
      public:
         /**
          * \brief Simple constructor
          *
          * \param spEqParams Shared equation parameters
          */
         ScalarEquation(SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~ScalarEquation();

         /**
          * @brief Set the shared pointer to the unknown field
          *
          * This is required because the field are not initialised at creation time
          *
          * \param spUnknown Shared pointer to the unknown of the equation
          */
         virtual void setUnknown(Code::SharedScalarVariable spUnknown);

         /**
          * @brief Get the unknown variable
          */
         const Code::ScalarVariable& unknown() const;

         /**
          * @brief Set the unknown variable
          */
         Code::ScalarVariable& rUnknown();

         /**
          * @brief Compute the nonlinear interaction term
          *
          * @param rNLComp Nonlinear term component
          */
         virtual void computeNonlinear(Code::PhysicalScalarType& rNLComp) = 0;

         /**
          * @brief Compute the linear term
          *
          * @param rRHS    RHS of timestepping equation
          */
         virtual void computeLinear(Code::SpectralScalarType& rRHS);

         /**
          * @brief Prepare the RHS for the timestep computation
          *
          * @param rhs    RHS of timestepping equation
          */
         virtual void prepareTimestep(const Code::SpectralScalarType& rhs) = 0;

         /**
          * @brief Transfer equation input to timestepper
          *
          * @param id      Component ID
          * @param storage Storage for the equation values
          * @param matIdx  Index of the given data
          * @param start   Start indx for the storage
          */
         virtual void copyTInput(FieldComponents::Spectral::Component id, DecoupledZMatrix& storage, const int matIdx, const int start);

         /**
          * @brief Transfer equation input to timestepper
          *
          * @param id      Component ID
          * @param storage Storage for the equation values
          * @param matIdx  Index of the given data
          * @param start   Start indx for the storage
          */
         virtual void copyTInput(FieldComponents::Spectral::Component id, MatrixZ& storage, const int matIdx, const int start);

         /**
          * @brief Transfer timestepper output to equation unknown
          *
          * @param id      Component ID
          * @param storage Storage for the equation values
          * @param matIdx  Index of the given data
          * @param start   Start indx for the storage
          */
         virtual void copyTOutput(FieldComponents::Spectral::Component id, const DecoupledZMatrix& storage, const int matIdx, const int start);

         /**
          * @brief Transfer timestepper output to equation unknown
          *
          * @param id      Component ID
          * @param storage Storage for the equation values
          * @param matIdx  Index of the given data
          * @param start   Start indx for the storage
          */
         virtual void copyTOutput(FieldComponents::Spectral::Component id, const MatrixZ& storage, const int matIdx, const int start);
         
      protected:
         /**
          * @brief Apply quasi-inverse to nonlinear values
          *
          * @param id      Component ID
          * @param storage Storage for the equation values
          * @param matIdx  Index of the given data
          * @param start   Start indx for the storage
          */
         virtual void applyQuasiInverse(FieldComponents::Spectral::Component id, DecoupledZMatrix& storage, const int matIdx, const int start);

         /**
          * @brief Apply quasi-inverse to nonlinear values
          *
          * @param id      Component ID
          * @param storage Storage for the equation values
          * @param matIdx  Index of the given data
          * @param start   Start indx for the storage
          */
         virtual void applyQuasiInverse(FieldComponents::Spectral::Component id, MatrixZ& storage, const int matIdx, const int start);

         /**
          * @brief Storage for the shared scalar variable
          */
         Code::SharedScalarVariable mspUnknown;

      private:
   };

   inline void ScalarEquation::setUnknown(Code::SharedScalarVariable spUnknown)
   {
      this->mspUnknown = spUnknown;
   }

   inline const Code::ScalarVariable& ScalarEquation::unknown() const
   {
      return *this->mspUnknown;
   }

   inline Code::ScalarVariable& ScalarEquation::rUnknown()
   {
      return *this->mspUnknown;
   }

   inline void ScalarEquation::computeLinear(Code::SpectralScalarType& rRHS)
   {
   }

   /// Typedef for a shared ScalarEquation
   typedef SharedPtrMacro<ScalarEquation> SharedScalarEquation;

}

#endif // SCALAREQUATION_HPP
