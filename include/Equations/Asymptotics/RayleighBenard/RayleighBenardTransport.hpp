/**
 * @file RayleighBenardTransport.hpp
 * @brief Implementation of the transport equation for the Rayleigh-Benard model 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 * @version 0.9.0
 * @date 2013-09-11
 */

#ifndef RAYLEIGHBENARDTRANSPORT_HPP
#define RAYLEIGHBENARDTRANSPORT_HPP

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
    * \brief Implementation of the transport equation for the Rayleigh-Benard model
    */
   class RayleighBenardTransport: public IScalarEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams  Shared equation parameters
          */
         RayleighBenardTransport(SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~RayleighBenardTransport();

         /**
          * @brief Compute the nonlinear interaction term
          *
          * @param rNLComp Nonlinear term component
          * @param id      ID of the component (allows for a more general implementation)
          */
         virtual void computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const;

         /**
          * @brief Set the equation matrices
          *
          * @param bcIds   List of boundary condition IDs
          * @param cbcIds  List of coupled boundary condition IDs
          */
         virtual void setSpectralMatrices(const BcEqMapType& bcIds, const std::map<PhysicalNames::Id, BcEqMapType>& cbcIds);
         
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

#endif // RAYLEIGHBENARDTRANSPORT_HPP
