/**
 * @file BoussinesqNoTiltedFPlane3DQGNoStreamfunction.hpp
 * @brief Implementation of the non orthogonal streamfunction computation for the Boussinesq tilted F-plane 3DQG model 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef BOUSSINESQNOTILTEDFPLANE3DQGNOSTREAMFUNTION_HPP
#define BOUSSINESQNOTILTEDFPLANE3DQGNOSTREAMFUNTION_HPP

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
    * @brief Implementation of the non orthogonal streamfunction computation for the Boussinesq tilted F-plane 3DQG model
    */
   class BoussinesqNoTiltedFPlane3DQGNoStreamfunction: public IScalarEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams Shared equation parameters
          */
         BoussinesqNoTiltedFPlane3DQGNoStreamfunction(SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~BoussinesqNoTiltedFPlane3DQGNoStreamfunction();

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

#endif // BOUSSINESQNOTILTEDFPLANE3DQGNOSTREAMFUNTION_HPP
