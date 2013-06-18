/** \file AnelasticFPlaneStreamfunction.hpp
 *  \brief Implementation of the streamfunction equation for the anelastic 3DQG f-plane model
 */

#ifndef ANELASTICFPLANESTREAMFUNCTION_HPP
#define ANELASTICFPLANESTREAMFUNCTION_HPP

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
    * \brief Implementation of the streamfunction equation for the anelastic 3DQG f-plane model
    */
   class AnelasticFPlaneStreamfunction: public IScalarEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams  Shared equation parameters
          */
         AnelasticFPlaneStreamfunction(SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~AnelasticFPlaneStreamfunction();

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
         virtual void setQuasiInverse(SparseMatrix &mat) const;

         /**
          * @brief Set the explicit linear matrix operator
          */
         virtual void setExplicitLinearBlock(DecoupledZSparse& mat, const SpectralFieldId fieldId, const MHDFloat k) const;

      private:
   };

   /**
    * @brief Get the quasi-inverse matrix operator
    *
    * @param eq      Equation to work on
    * @param mat     Storage for output matrix
    */
   void quasiInverseBlock(const AnelasticFPlaneStreamfunction& eq, SparseMatrix& mat);

   /**
    * @brief Get the time matrix block
    *
    * @param eq      Equation to work on
    * @param mat     Storage for output matrix
    * @param k       Wave number k
    */
   void timeBlock(const AnelasticFPlaneStreamfunction& eq, DecoupledZSparse& mat, const MHDFloat k);

   /**
    * @brief Get the linear matrix block on given field
    *
    * @param eq      Equation to work on
    * @param mat     Storage for output matrix
    * @param fieldId Physical ID of the field
    * @param k       Wave number k
    */
   void linearBlock(const AnelasticFPlaneStreamfunction& eq, DecoupledZSparse& mat, const SpectralFieldId fieldId, const MHDFloat k);

   /**
    * @brief Get the boundary condition matrix block on given field
    *
    * @param eq      Equation to work on
    * @param mat     Storage for output matrix
    * @param fieldId Physical ID of the field
    * @param k       Wave number k
    */
   void boundaryBlock(const AnelasticFPlaneStreamfunction& eq, DecoupledZSparse& mat, const SpectralFieldId fieldId, const MHDFloat k);

}
}

#endif // ANELASTICFPLANESTREAMFUNCTION_HPP
