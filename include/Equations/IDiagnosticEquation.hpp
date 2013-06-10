/** \file IDiagnosticEquation.hpp
 *  \brief Base building block for the implementation of a diagnostic equation
 */

#ifndef IDIAGNOSTICEQUATION_HPP
#define IDIAGNOSTICEQUATION_HPP

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
#include "Enums/Dimensions.hpp"
#include "Enums/FieldIds.hpp"
#include "SpectralOperators/BoundaryConditions.hpp"
#include "Equations/EquationParameters.hpp"
#include "Equations/CouplingInformation.hpp"
#include "Equations/EquationData.hpp"
#include "TypeSelectors/VariableSelector.hpp"
#include "Variables/VariableRequirement.hpp"
#include "Simulation/SimulationBoundary.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   /**
    * @brief Base building block for the implementation of a diagnostic equation
    */
   class IDiagnosticEquation : public EquationData
   {
      public:
         /**
          * @brief Simple constructor
          */
         explicit IDiagnosticEquation(SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~IDiagnosticEquation();

         /**
          * @brief Initialise the equation
          */
         virtual void init();

         /**
          * @brief Compute the explicit linear terms
          *
          * @param compId     Equation field component ID
          * @param eqField    Equation field values
          * @param eqStart    Start index for the equation field
          * @param fieldId    Physical field ID
          * @param linField   Explicit linear field values
          * @param linStart   Start index for the linear field
          * @param matIdx     System index
          */
         virtual void computeLinear(FieldComponents::Spectral::Id compId, DecoupledZMatrix& eqField, const int eqStart, SpectralFieldId fieldId, const DecoupledZMatrix& linField, const int linStart, const int matIdx) const;
         virtual void computeLinear(FieldComponents::Spectral::Id compId, DecoupledZMatrix& eqField, const int eqStart, SpectralFieldId fieldId, const MatrixZ& linField, const int linStart, const int matIdx) const;

         /**
          * @brief Compute the explicit linear terms
          *
          * @param compId     Equation field component ID
          * @param eqField    Equation field values
          * @param eqStart    Start index for the equation field
          * @param linField   Explicit linear field values
          * @param linStart   Start index for the linear field
          * @param matIdx     System index
          */
         virtual void computeLinear(FieldComponents::Spectral::Id compId, MatrixZ& eqField, const int eqStart, SpectralFieldId fieldId, const DecoupledZMatrix& linField, const int linStart, const int matIdx) const;
         virtual void computeLinear(FieldComponents::Spectral::Id compId, MatrixZ& eqField, const int eqStart, SpectralFieldId fieldId, const MatrixZ& linField, const int linStart, const int matIdx) const;

         /**
          * @brief Apply quasi-inverse to the nonlinear term
          *
          * @param id      Component ID
          * @param storage Storage for the equation values
          * @param matIdx  Index of the given data
          * @param start   Start index for the storage
          */
         virtual void applyNLQuasiInverse(FieldComponents::Spectral::Id id, DecoupledZMatrix& storage, const int matIdx, const int start);

         /**
          * @brief Apply quasi-inverse to the nonlinear term
          *
          * @param id      Component ID
          * @param storage Storage for the equation values
          * @param matIdx  Index of the given data
          * @param start   Start index for the storage
          */
         virtual void applyNLQuasiInverse(FieldComponents::Spectral::Id id, MatrixZ& storage, const int matIdx, const int start);
         
      protected:
         /**
          * @brief Set the equation variable requirements
          */
         virtual void setRequirements() = 0;

         /**
          * @brief Set the equation coupling information
          */
         virtual void setCoupling() = 0;

      private:
   };

   /// Typedef for a smart IDiagnosticEquation
   typedef SharedPtrMacro<IDiagnosticEquation> SharedIDiagnosticEquation;
}
}

#endif // IDIAGNOSTICEQUATION_HPP
