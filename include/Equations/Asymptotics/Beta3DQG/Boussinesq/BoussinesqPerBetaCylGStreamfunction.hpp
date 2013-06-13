/** \file BoussinesqPerBetaCylGStreamfunction.hpp
 *  \brief Implementation of the streamfunction equation for the Boussinesq beta model with cylindrical gravity with periodic radius
 */

#ifndef BOUSSINESQPERBETACYLGSTREAMFUNCTION_HPP
#define BOUSSINESQPERBETACYLGSTREAMFUNCTION_HPP

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
#include "Equations/Asymptotics/Beta3DQG/Boussinesq/IBoussinesqPerBetaCylGScalarEquation.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   /**
    * \brief Implementation of the streamfunction equation for the Boussinesq beta model with cylindrical gravity with periodic radius
    */
   class BoussinesqPerBetaCylGStreamfunction: public IBoussinesqPerBetaCylGScalarEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams  Shared equation parameters
          */
         BoussinesqPerBetaCylGStreamfunction(SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~BoussinesqPerBetaCylGStreamfunction();

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
         
      protected:
         /**
          * @brief Set variable requirements
          */
         virtual void setRequirements();

         /**
          * @brief Set the equation coupling information
          */
         virtual void setCoupling();

      private:
   };

   /**
    * @brief Get the time matrix block for an equation
    *
    * @param mat     Storage for output matrix
    * @param nX      Matrix size in X
    * @param nZ      Matrix size in Z
    * @param k       Wave number k
    */
   void timeBlock(const BoussinesqBetaCylGStreamfunction& eq, DecoupledZSparse& mat, const int nX, const int nZ, const MHDFloat k);

   /**
    * @brief Get the linear matrix block for an equation on given field
    *
    * @param mat     Storage for output matrix
    * @param fieldId Physical ID of the field
    * @param nX      Matrix size in X
    * @param nZ      Matrix size in Z
    * @param k       Wave number k
    */
   void linearBlock(const BoussinesqBetaCylGStreamfunction& eq, DecoupledZSparse& mat, const SpectralFieldId fieldId, const int nX, const int nZ, const MHDFloat k);

   /**
    * @brief Get the boundary condition matrix block for an equation on given field
    *
    * @param mat     Storage for output matrix
    * @param fieldId Physical ID of the field
    * @param nX      Matrix size in X
    * @param nZ      Matrix size in Z
    * @param k       Wave number k
    */
   void boundaryBlock(const BoussinesqBetaCylGStreamfunction& eq, DecoupledZSparse& mat, const SpectralFieldId fieldId, const int nX, const int nZ, const MHDFloat k);

}
}

#endif // BOUSSINESQPERBETACYLGSTREAMFUNCTION_HPP
