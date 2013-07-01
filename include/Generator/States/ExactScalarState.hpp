/** \file ExactScalarState.hpp
 *  \brief Implementation of the equation to generate exact scalar states
 */

#ifndef EXACTSCALARSTATE_HPP
#define EXACTSCALARSTATE_HPP

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
    * @brief Implementation of the equation to generate exact scalar state
    */
   class ExactScalarState: public IScalarEquation
   {
      public:
         /**
          * @brief Enums for the avaialable exact states
          */
         enum StateTypeId {
            CONSTANT,
            SINESINE,
            SINECOSINE
         };

         /**
          * @brief Simple constructor
          *
          * @param spEqParams Shared equation parameters
          * @paarm name       Name of the field
          */
         ExactScalarState(SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~ExactScalarState();

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
         virtual DecoupledZSparse operatorRow(const OperatorRowId opId, FieldComponents::Spectral::Id comp, const int matIdx) const;

         /**
          * @brief Initialise spectral equation matrices
          *
          * @param spBcIds   List of boundary condition IDs
          */
         virtual void initSpectralMatrices(const SharedSimulationBoundary spBcIds);

         /**
          * @brief Set the unknown name and requirements
          */
         void setIdentity(const PhysicalNames::Id name);

         /**
          * @brief Set the state type id
          */
         void setStateType(const ExactScalarState::StateTypeId id);

         /**
          * @brief Set the options for the sine state
          *
          * @param aX   Amplitude of the sine in X direction
          * @param kX   Wave number of the sine in X direction
          * @param aZ   Amplitude of the sine in Z direction
          * @param kZ   Wave number of the sine in Z direction
          */
         void setSineOptions(const MHDFloat aX, const MHDFloat kX, const MHDFloat aZ, const MHDFloat kZ);

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
         virtual void setExplicitLinearBlock(FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const MHDFloat k) const;

      private:
         /**
          * @brief Type of the state to generate
          */
         StateTypeId mTypeId;

         /**
          * @brief Amplitude of the sine state
          */
         Array mSineA;

         /**
          * @brief Wave number of the sine state
          */
         Array mSineN;
   };

   /// Typedef for a shared ExactScalarState
   typedef SharedPtrMacro<ExactScalarState> SharedExactScalarState;

   /**
    * @brief Get the quasi-inverse matrix operator
    *
    * @param eq      Equation to work on
    * @param mat     Storage for output matrix
    */
   void quasiInverseBlock(const ExactScalarState& eq, FieldComponents::Spectral::Id compId, SparseMatrix& mat);

   /**
    * @brief Get the linear matrix block for an equation on given field
    *
    * @param mat     Storage for output matrix
    * @param fieldId Physical ID of the field
    * @param k       Wave number k
    */
   void linearBlock(const ExactScalarState& eq, FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const MHDFloat k);

   /**
    * @brief Get the boundary condition matrix block for an equation on given field
    *
    * @param mat     Storage for output matrix
    * @param fieldId Physical ID of the field
    * @param k       Wave number k
    */
   void boundaryBlock(const ExactScalarState& eq, FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const MHDFloat k);

}
}

#endif // EXACTSCALARSTATE_HPP
