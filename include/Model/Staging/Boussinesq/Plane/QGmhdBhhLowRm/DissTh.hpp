/**
 * @file DissTh.hpp
 * @brief Implementation of the thermal dissipation equation for the Boussinesq F-plane QG model with horizontal helicoidal magnetic field applied
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 * @modified by Stefano Maffei \<maffei.ste@gmail.com\>
 */

#ifndef QUICC_EQUATIONS_BOUSSINESQ_PLANE_QGMHDBHHLOWRM_DISSTH_HPP
#define QUICC_EQUATIONS_BOUSSINESQ_PLANE_QGMHDBHHLOWRM_DISSTH_HPP

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

namespace QuICC {

namespace Equations {

namespace Boussinesq {

namespace Plane {

namespace QGmhdBhhLowRm {
   /**
    * @brief Implementation of the thermal dissipation equation for the Boussinesq F-plane QG model
    */
   class DissTh: public IScalarEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams Shared equation parameters
          */
         DissTh(SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~DissTh();
         
         /**
          * @brief Compute the nonlinear interaction term
          *
          * @param rNLComp Nonlinear term component
          * @param id      ID of the component (allows for a more general implementation)
          */
         virtual void computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const;
         
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

}
}
}
}
}
#endif // QUICC_EQUATIONS_BOUSSINESQ_PLANE_QGMHDBHHLOWRM_DISSTH_HPP
