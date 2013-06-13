/** \file IBoussinesqBetaCylGScalarEquation.hpp
 *  \brief Implementation of the general scalar equation for the Boussinesq beta model with cylindrical gravity
 */

#ifndef IBOUSSINESQBETACYLGSCALAREQUATION_HPP
#define IBOUSSINESQBETACYLGSCALAREQUATION_HPP

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
#include "Simulation/SimulationBoundary.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   /**
    * @brief Implementation of the general scalar equation for the Boussinesq beta model with cylindrical gravity
    */
   class IBoussinesqBetaCylGScalarEquation: public IScalarEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams  Shared equation parameters
          */
         IBoussinesqBetaCylGScalarEquation(SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~IBoussinesqBetaCylGScalarEquation();

         /**
          * @brief Initialise spectral equation matrices
          *
          * @param spBcs   List of boundary condition IDs
          */
         virtual void initSpectralMatrices(const SharedSimulationBoundary bcs);
         
      protected:
         /**
          * @brief Storage for the shared boundary condition list
          */
         SharedSimulationBoundary mspBcIds;

      private:
   };

}
}

#endif // IBOUSSINESQBETACYLGSCALAREQUATION_HPP
