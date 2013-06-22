/** \file BoussinesqPerPerBetaCylGVertical.hpp
 *  \brief Implementation of the vertical velocity equation for the Boussinesq beta model with cylindrical gravity with periodic radius
 */

#ifndef BOUSSINESQPERBETACYLGVERTICAL_HPP
#define BOUSSINESQPERBETACYLGVERTICAL_HPP

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
    * \brief Implementation of the vertical velocity equation for the Boussinesq beta model with cylindrical gravity with periodic radius
    */
   class BoussinesqPerPerBetaCylGVertical: public IScalarEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams  Shared equation parameters
          */
         BoussinesqPerPerBetaCylGVertical(SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~BoussinesqPerPerBetaCylGVertical();

         /**
          * @brief Compute the nonlinear interaction term
          *
          * @param rNLComp Nonlinear term component
          * @param id      ID of the component (allows for a more general implementation)
          */
         virtual void computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const;

         /**
          * @brief Generic operator row dispatcher
          */
         virtual DecoupledZSparse operatorRow(const OperatorRowId opId, FieldComponents::Spectral::Id comp, const int matIdx) const;

         /**
          * @brief Initialise the spectral equation matrices
          *
          * @param spBcIds   List of boundary condition IDs
          */
         virtual void initSpectralMatrices(const SharedSimulationBoundary spBcIds);
         
      protected:
         /**
          * @brief Set variable requirements
          */
         virtual void setRequirements();

         /**
          * @brief Set the equation coupling information
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
   };

   /**
    * @brief Get the quasi-inverse matrix operator
    *
    * @param eq      Equation to work on
    * @param mat     Storage for output matrix
    */
   void quasiInverseBlock(const BoussinesqPerBetaCylGStreamfunction& eq, SparseMatrix& mat);

   /**
    * @brief Get the time matrix block for an equation
    *
    * @param mat     Storage for output matrix
    * @param nX      Matrix size in X
    * @param nZ      Matrix size in Z
    * @param k       Wave number k
    */
   void timeBlock(const BoussinesqPerBetaCylGStreamfunction& eq, DecoupledZSparse& mat, const int nX, const int nZ, const MHDFloat k);

   /**
    * @brief Get the linear matrix block for an equation on given field
    *
    * @param mat     Storage for output matrix
    * @param fieldId Physical ID of the field
    * @param nX      Matrix size in X
    * @param nZ      Matrix size in Z
    * @param k       Wave number k
    */
   void linearBlock(const BoussinesqPerBetaCylGStreamfunction& eq, DecoupledZSparse& mat, const SpectralFieldId fieldId, const int nX, const int nZ, const MHDFloat k);

   /**
    * @brief Get the boundary condition matrix block for an equation on given field
    *
    * @param mat     Storage for output matrix
    * @param fieldId Physical ID of the field
    * @param nX      Matrix size in X
    * @param nZ      Matrix size in Z
    * @param k       Wave number k
    */
   void boundaryBlock(const BoussinesqPerBetaCylGStreamfunction& eq, DecoupledZSparse& mat, const SpectralFieldId fieldId, const int nX, const int nZ, const MHDFloat k);

}
}

#endif // BOUSSINESQPERBETACYLGVERTICAL_HPP
