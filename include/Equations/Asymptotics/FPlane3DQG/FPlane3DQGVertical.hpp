/** \file FPlane3DQGVertical.hpp
 *  \brief Implementation of the vertical velocity equation for the 3DQG f-plane model
 */

#ifndef FPLANE3DQGVERTICAL_HPP
#define FPLANE3DQGVERTICAL_HPP

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
#include "Equations/IScalarPEquation.hpp"
#include "TypeSelectors/ScalarSelector.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   /**
    * \brief Implementation of the vertical velocity equation for the 3DQG f-plane model
    */
   class FPlane3DQGVertical: public IScalarPEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams  Shared equation parameters
          */
         FPlane3DQGVertical(SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~FPlane3DQGVertical();

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

#endif // FPLANE3DQGVERTICAL_HPP
