/** \file IScalarPEquation.hpp
 *  \brief Base for the implementation of a scalar prognostic equation
 */

#ifndef ISCALARPEQUATION_HPP
#define ISCALARPEQUATION_HPP

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
#include "Equations/IPrognosticEquation.hpp"
#include "TypeSelectors/ScalarSelector.hpp"
#include "TypeSelectors/VariableSelector.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   /**
    * \brief Base for the implementation of a scalar prognostic equation
    */
   class IScalarPEquation: public IPrognosticEquation
   {
      public:
         /**
          * \brief Simple constructor
          *
          * \param spEqParams Shared equation parameters
          */
         explicit IScalarPEquation(SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~IScalarPEquation();

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

         /**
          * @brief Compute the nonlinear interaction term
          *
          * @param rNLComp Nonlinear term component
          */
         virtual void computeNonlinear(Datatypes::PhysicalScalarType& rNLComp) const = 0;

         /**
          * @brief Prepare the RHS for the timestep computation
          *
          * @param rhs    RHS of timestepping equation
          */
         virtual void prepareTimestep(const Datatypes::SpectralScalarType& rhs);

         /**
          * @brief Basic implementatin of the method to transfer equation input to timestepper
          *
          * @param id         Component ID
          * @param storage    Storage for the equation values
          * @param matIdx     Index of the given data
          * @param start      Start index for the storage
          */
         virtual void timestepInput(FieldComponents::Spectral::Id id, DecoupledZMatrix& storage, const int matIdx, const int start);

         /**
          * @brief Basic implementation of the method to transfer equation input to timestepper
          *
          * @param id         Component ID
          * @param storage    Storage for the equation values
          * @param matIdx     Index of the given data
          * @param start      Start index for the storage
          */
         virtual void timestepInput(FieldComponents::Spectral::Id id, MatrixZ& storage, const int matIdx, const int start);

         /**
          * @brief Basic implementation of the method to transfer timestepper output to equation unknown
          *
          * @param id      Component ID
          * @param storage Storage for the equation values
          * @param matIdx  Index of the given data
          * @param start   Start index for the storage
          */
         virtual void timestepOutput(FieldComponents::Spectral::Id id, const DecoupledZMatrix& storage, const int matIdx, const int start);

         /**
          * @brief Basic implementation of the method to transfer timestepper output to equation unknown
          *
          * @param id      Component ID
          * @param storage Storage for the equation values
          * @param matIdx  Index of the given data
          * @param start   Start index for the storage
          */
         virtual void timestepOutput(FieldComponents::Spectral::Id id, const MatrixZ& storage, const int matIdx, const int start);

         /**
          * @brief Transfer equation input to timestepper
          *
          * @param id      Component ID
          * @param storage Storage for the equation values
          * @param matIdx  Index of the given data
          * @param start   Start index for the storage
          */
         virtual void copyTInput(FieldComponents::Spectral::Id id, DecoupledZMatrix& storage, const int matIdx, const int start);

         /**
          * @brief Transfer equation input to timestepper
          *
          * @param id      Component ID
          * @param storage Storage for the equation values
          * @param matIdx  Index of the given data
          * @param start   Start index for the storage
          */
         virtual void copyTInput(FieldComponents::Spectral::Id id, MatrixZ& storage, const int matIdx, const int start);

         /**
          * @brief Transfer timestepper output to equation unknown
          *
          * @param id      Component ID
          * @param storage Storage for the equation values
          * @param matIdx  Index of the given data
          * @param start   Start index for the storage
          */
         virtual void copyTOutput(FieldComponents::Spectral::Id id, const DecoupledZMatrix& storage, const int matIdx, const int start);

         /**
          * @brief Transfer timestepper output to equation unknown
          *
          * @param id      Component ID
          * @param storage Storage for the equation values
          * @param matIdx  Index of the given data
          * @param start   Start index for the storage
          */
         virtual void copyTOutput(FieldComponents::Spectral::Id id, const MatrixZ& storage, const int matIdx, const int start);
         
      protected:

      private:
         /**
          * @brief Storage for the shared scalar variable
          */
         Datatypes::SharedScalarVariableType mspUnknown;
   };

   /// Typedef for a shared IScalarPEquation
   typedef SharedPtrMacro<IScalarPEquation> SharedIScalarPEquation;

}
}

#endif // ISCALARPEQUATION_HPP
