/** \file IDiagnosticEquation.hpp
 *  \brief Base building block for the implementation of a diagnostic equation
 */

#ifndef IDIAGNOSTICEQUATION_HPP
#define IDIAGNOSTICEQUATION_HPP

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
#include "Enums/Dimensions.hpp"
#include "Enums/FieldIds.hpp"
#include "SpectralOperators/BoundaryConditions.hpp"
#include "Equations/EquationParameters.hpp"
#include "Equations/CouplingInformation.hpp"
#include "Equations/IEquation.hpp"
#include "TypeSelectors/VariableSelector.hpp"
#include "Variables/VariableRequirement.hpp"
#include "Simulation/SimulationBoundary.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   /**
    * @brief Base building block for the implementation of a diagnostic equation
    */
   class IDiagnosticEquation : public IEquation
   {
      public:
         /**
          * @brief Simple constructor
          */
         explicit IDiagnosticEquation(SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~IDiagnosticEquation();
        
      protected:

      private:
   };

   /// Typedef for a smart IDiagnosticEquation
   typedef SharedPtrMacro<IDiagnosticEquation> SharedIDiagnosticEquation;
}
}

#endif // IDIAGNOSTICEQUATION_HPP
