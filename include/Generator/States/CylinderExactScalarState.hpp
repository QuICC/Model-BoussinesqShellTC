/**
 * @file CylinderExactScalarState.hpp
 * @brief Implementation of the equation to generate exact scalar states in a cylinder 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef CYLINDEREXACTSCALARSTATE_HPP
#define CYLINDEREXACTSCALARSTATE_HPP

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
    * @brief Implementation of the equation to generate exact scalar state in a cylinder
    */
   class CylinderExactScalarState: public IScalarEquation
   {
      public:
         /**
          * @brief Enums for the avaialable exact states
          */
         enum StateTypeId {
            CONSTANT,
         };

         /**
          * @brief Simple constructor
          *
          * @param spEqParams Shared equation parameters
          * @paarm name       Name of the field
          */
         CylinderExactScalarState(SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~CylinderExactScalarState();

         /**
          * @brief Compute the nonlinear interaction term
          *
          * @param rNLComp Nonlinear term component
          * @param compId  ID of the component (allows for a more general implementation)
          */
         virtual void computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id compId) const;

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
          * @brief Generic operator row dispatcher
          */
         virtual DecoupledZSparse operatorRow(const OperatorRowId opId, FieldComponents::Spectral::Id comp, const int matIdx, const bool hasBoundary) const;

         /**
          * @brief Set the unknown name and requirements
          */
         void setIdentity(const PhysicalNames::Id name);

         /**
          * @brief Set the state type id
          */
         void setStateType(const CylinderExactScalarState::StateTypeId id);

      protected:
         /**
          * @brief Set variable requirements
          */
         virtual void setRequirements();

         /**
          * @brief Set coupling information
          */
         virtual void setCoupling();

         /**
          * @brief Set the quasi inverse matrix operator
          */
         virtual void setQuasiInverse(FieldComponents::Spectral::Id compId, SparseMatrix &mat) const;

         /**
          * @brief Set the explicit linear matrix operator
          */
         virtual void setExplicitLinearBlock(FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const std::vector<MHDFloat>& eigs) const;

      private:
         /**
          * @brief Type of the state to generate
          */
         StateTypeId mTypeId;
   };

   /// Typedef for a shared CylinderExactScalarState
   typedef SharedPtrMacro<CylinderExactScalarState> SharedCylinderExactScalarState;

   /**
    * @brief Get the quasi-inverse matrix operator
    *
    * @param eq      Equation to work on
    * @param mat     Storage for output matrix
    */
   void quasiInverseBlock(const CylinderExactScalarState& eq, FieldComponents::Spectral::Id compId, SparseMatrix& mat);

   /**
    * @brief Get the linear matrix block for an equation on given field
    *
    * @param mat     Storage for output matrix
    * @param fieldId Physical ID of the field
    * @parm eigs     Wave number k
    */
   void linearBlock(const CylinderExactScalarState& eq, FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const std::vector<MHDFloat>& eigs, const bool hasBoundary);

   /**
    * @brief Get the boundary condition matrix block for an equation on given field
    *
    * @param eq      Equation to work on
    * @param compId  Component ID
    * @param fieldId Physical ID of the field
    * @parm eigs     Wave number k
    */
   void boundaryBlock(CylinderExactScalarState& eq, FieldComponents::Spectral::Id compId, const SpectralFieldId fieldId, const std::vector<MHDFloat>& eigs, std::vector<MHDFloat>& coeffs, std::vector<Boundary::BCIndex>& bcIdx);

}
}

#endif // CYLINDEREXACTSCALARSTATE_HPP