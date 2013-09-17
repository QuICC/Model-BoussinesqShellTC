/**
 * @file IEquation.hpp
 * @brief Base building block for the implementation of an equation 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef IEQUATION_HPP
#define IEQUATION_HPP

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
    * @brief Base building block for the implementation of an equation
    */
   class IEquation : public EquationData
   {
      public:
         /**
          * @brief Enum for the different types of operator rows
          */
         enum OperatorRowId {
            TIMEROW = 0,
            LINEARROW,
            BOUNDARYROW};

         /**
          * @brief Simple constructor
          */
         explicit IEquation(SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~IEquation();

         /**
          * @brief Access the shared resolution
          */
         virtual SharedResolution spRes() const = 0;

         /**
          * @brief Initialise the equation
          */
         virtual void init();

         /**
          * @brief Generic operator row dispatcher
          *
          * It has only a dummy implementation and should never get called
          */
         virtual DecoupledZSparse  operatorRow(const OperatorRowId opId, FieldComponents::Spectral::Id comp, const int matIdx) const; //= 0;

         /**
          * @brief Compute the nonlinear interaction term
          *
          * @param rNLComp Nonlinear term component
          * @param id      ID of the component (allows for a more general implementation)
          */
         virtual void computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const;

         /**
          * @brief Compute the source term
          *
          * @param compId  ID of the spectral component
          * @param i       Fastest index
          * @param j       Second index
          * @param k       Slowest index
          */
         virtual MHDComplex sourceTerm(FieldComponents::Spectral::Id compId, const int i, const int j, const int k) const;

         /**
          * @brief Initialise the spectral equation matrices
          *
          * @param spBcIds   List of boundary condition IDs
          */
         virtual void initSpectralMatrices(const SharedSimulationBoundary spBcIds) = 0;
         
      protected:
         /**
          * @brief Set the equation variable requirements
          */
         virtual void setRequirements() = 0;

         /**
          * @brief Set the equation coupling information
          */
         virtual void setCoupling() = 0;

         /**
          * @brief Initialise the spectral equation matrices with a single "eigen" direction
          *
          * @param spBcIds List of boundary condition IDs
          * @param compId  Spectral component
          */
         void initSpectralMatrices1DEigen(const SharedSimulationBoundary spBcIds, FieldComponents::Spectral::Id compId);

         /**
          * @brief Initialise the spectral equation matrices with two "eigen" direction
          *
          * @param spBcIds List of boundary condition IDs
          * @param compId  Spectral component
          */
         void initSpectralMatrices2DEigen(const SharedSimulationBoundary spBcIds, FieldComponents::Spectral::Id compId);

      private:
         /**
          * @brief Set the quasi inverse matrix operator
          *
          * It has only a dummy implementation and should never get called!
          */
         virtual void setQuasiInverse(FieldComponents::Spectral::Id compId, SparseMatrix &mat) const; //= 0;

         /**
          * @brief Set the explicit linear matrix operator with a single eigen direction
          */
         virtual void setExplicitLinearBlock(FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const MHDFloat k) const; //= 0;

         /**
          * @brief Set the explicit linear matrix operator with two eigen direction
          */
         virtual void setExplicitLinearBlock(FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const MHDFloat k2D, const MHDFloat k3D) const; //= 0;

   };

   /// Typedef for a smart IEquation
   typedef SharedPtrMacro<IEquation> SharedIEquation;

   /**
    * @brief Apply quasi-inverse to the nonlinear term
    *
    * @param compId  Component ID
    * @param storage Storage for the equation values
    * @param matIdx  Index of the given data
    * @param start   Start index for the storage
    */
   void applyQuasiInverse(const IEquation& eq, FieldComponents::Spectral::Id compId, DecoupledZMatrix& storage, const int matIdx, const int start);
   void applyQuasiInverse(const IEquation& eq, FieldComponents::Spectral::Id compId, MatrixZ& storage, const int matIdx, const int start);

   /**
    * @brief Compute and add the explicit linear terms
    *
    * @param compId     Equation field component ID
    * @param eqField    Equation field values
    * @param eqStart    Start index for the equation field
    * @param fieldId    Physical field ID
    * @param explicitField   Explicit linear field values
    * @param matIdx     System index
    */
   void addExplicitLinear(const IEquation& eq, FieldComponents::Spectral::Id compId, DecoupledZMatrix& eqField, const int eqStart, SpectralFieldId fieldId, const Datatypes::SpectralScalarType& explicitField, const int matIdx);
   void addExplicitLinear(const IEquation& eq, FieldComponents::Spectral::Id compId, MatrixZ& eqField, const int eqStart, SpectralFieldId fieldId, const Datatypes::SpectralScalarType& explicitField, const int matIdx);

   /**
    * @brief Dummy implementation: This should never be called!
    */
   void quasiInverseBlock(const IEquation& eq, SparseMatrix& mat); //= 0;

   /**
    * @brief Dummy implementation: This should never be called!
    */
   void timeBlock(const IEquation& eq, FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const MHDFloat k); //= 0;

   /**
    * @brief Dummy implementation: This should never be called!
    */
   void linearBlock(const IEquation& eq, FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const MHDFloat k); //= 0;

   /**
    * @brief Dummy implementation. This should never get called!
    */
   void boundaryBlock(const IEquation& eq, FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const MHDFloat k); //= 0;
   
}
}

#endif // IEQUATION_HPP
