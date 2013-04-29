/** \file BoussinesqBetaSphGStreamfunction.hpp
 *  \brief Implementation of the streamfunction equation for the Boussinesq beta model with spherical gravity
 */

#ifndef BOUSSINESQBETASPHGSTREAMFUNCTION_HPP
#define BOUSSINESQBETASPHGSTREAMFUNCTION_HPP

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
#include "Equations/Asymptotics/Beta3DQG/IBoussinesqBetaSphGScalarEquation.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   /**
    * \brief Implementation of the streamfunction equation for the Boussinesq beta model with spherical gravity
    */
   class BoussinesqBetaSphGStreamfunction: public IBoussinesqBetaSphGScalarEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams  Shared equation parameters
          */
         BoussinesqBetaSphGStreamfunction(SharedIEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~BoussinesqBetaSphGStreamfunction();

         /**
          * @brief Compute the nonlinear interaction term
          *
          * @param rNLComp Nonlinear term component
          */
         virtual void computeNonlinear(Datatypes::PhysicalScalarType& rNLComp) const;

         /**
          * @brief Transfer timestepper output to unknown and update vorticity
          *
          * @param id      Component ID
          * @param storage Storage for the equation values
          * @param matIdx  Index of the given data
          * @param start   Start indx for the storage
          */
         virtual void timestepOutput(FieldComponents::Spectral::Id id, const DecoupledZMatrix& storage, const int matIdx, const int start);

         /**
          * @brief Transfer timestepper output to unknown and update vorticity
          *
          * @param id      Component ID
          * @param storage Storage for the equation values
          * @param matIdx  Index of the given data
          * @param start   Start indx for the storage
          */
         virtual void timestepOutput(FieldComponents::Spectral::Id id, const MatrixZ& storage, const int matIdx, const int start);
         
      protected:
         /**
          * @brief Set variable requirements
          */
         virtual void setRequirements();

      private:
   };

}
}

#endif // BOUSSINESQBETASPHGSTREAMFUNCTION_HPP
