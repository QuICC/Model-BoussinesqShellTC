/** \file RotatingRBStreamfunction.hpp
 *  \brief Implementation of the streamfunction equation for the rotating Rayleigh-Benard model
 */

#ifndef ROTATINGRBSTREAMFUNCTION_HPP
#define ROTATINGRBSTREAMFUNCTION_HPP

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
    * \brief Implementation of the streamfunction equation for the rotating Rayleigh-Benard model
    */
   class RotatingRBStreamfunction: public IScalarEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams  Shared equation parameters
          */
         RotatingRBStreamfunction(SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~RotatingRBStreamfunction();

         /**
          * @brief Compute the nonlinear interaction term
          *
          * @param rNLComp Nonlinear term component
          * @param id      ID of the component (allows for a more general implementation)
          */
         virtual void computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const;

         /**
          * @brief Set equation matrices
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

#endif // ROTATINGRBSTREAMFUNCTION_HPP
