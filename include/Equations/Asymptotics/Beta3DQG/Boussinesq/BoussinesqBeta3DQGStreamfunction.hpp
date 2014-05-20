/**
 * @file BoussinesqBeta3DQGStreamfunction.hpp
 * @brief Implementation of the streamfunction equation for the Boussinesq Beta 3DQG model
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef BOUSSINESQBETA3DQGSTREAMFUNCTION_HPP
#define BOUSSINESQBETA3DQGSTREAMFUNCTION_HPP

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
    * @brief Implementation of the streamfunction equation for the Boussinesq Beta 3DQG model
    */
   class BoussinesqBeta3DQGStreamfunction: public IScalarEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param pyName     Python script name
          * @param spEqParams Shared equation parameters
          */
         BoussinesqBeta3DQGStreamfunction(const std::string& pyName, SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~BoussinesqBeta3DQGStreamfunction();

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

#endif // BOUSSINESQBETA3DQGSTREAMFUNCTION_HPP
