/** \file Beta3DQGTransport.hpp
 *  \brief Implementation of the transport equation for the 3DQG beta model
 */

#ifndef BETA3DQGTRANSPORT_HPP
#define BETA3DQGTRANSPORT_HPP

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
#include "Equations/IScalarEquation.hpp"
#include "TypeSelectors/ScalarSelector.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   /**
    * \brief Implementation of the transport equation for the 3DQG beta model
    */
   class Beta3DQGTransport: public IScalarEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams  Shared equation parameters
          */
         Beta3DQGTransport(SharedIEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~Beta3DQGTransport();

         /**
          * @brief Compute the nonlinear interaction term
          *
          * @param rNLComp Nonlinear term component
          */
         virtual void computeNonlinear(Datatypes::PhysicalScalarType& rNLComp) const;

         /**
          * @brief Compute the linear term
          *
          * @param rRHS    RHS of timestepping equation
          */
         virtual void computeLinear(Datatypes::SpectralScalarType& rRHS) const;

         /**
          * @brief Build Full block row for linear operators
          */
         virtual DecoupledZSparse linearRow(FieldComponents::Spectral::Id comp, const int matIdx) const;

         /**
          * @brief Build Full block row for time operators
          */
         virtual DecoupledZSparse timeRow(FieldComponents::Spectral::Id comp, const int matIdx) const;

         /**
          * @brief Build Full block row for time operators
          */
         virtual DecoupledZSparse boundaryRow(FieldComponents::Spectral::Id comp, const int matIdx) const;

         /**
          * @brief Initialise spectral equation matrices
          *
          * @param spBcs   List of boundary condition IDs
          */
         virtual void initSpectralMatrices(const SharedSimulationBoundary bcs);
         
      protected:
         /**
          * @brief Set variable requirements
          */
         virtual void setRequirements();

         /**
          * @brief Set coupling information
          */
         virtual void setCoupling();

      private:
   };

}
}

#endif // BETA3DQGTRANSPORT_HPP
