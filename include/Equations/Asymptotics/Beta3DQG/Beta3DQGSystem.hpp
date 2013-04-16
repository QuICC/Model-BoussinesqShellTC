/** \file Beta3DQGSystem.hpp
 *  \brief Implementation of the system of equations equation for the 3DQG beta model
 */

#ifndef BETA3DQGSYSTEM_HPP
#define BETA3DQGSYSTEM_HPP

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

namespace GeoMHDiSCC {

namespace Equations {

   /**
    *  @brief Implementation of the system of equations equation for the 3DQG beta model
    */
   class Beta3DQGSystem
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams  Shared equation parameters
          */
         Beta3DQGSystem(SharedIEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~Beta3DQGSystem();

         const std::vector<PhysicalNames::Id> fields(const PhysicalNames::Id eqId);
         void linear(DecoupledZSparse& op, const MHDFloat k, const PhysicalNames::Id eqId, const PhysicalNames::Id fieldId);
         void time(DecoupledZSparse& op, const PhysicalNames::Id eqId);
         void boundary(DecoupledZSparse& op, const PhysicalNames::Id eqId, const PhysicalNames::Id fieldId);
         void quasiInverse(DecoupledZSparse& op, const PhysicalNames::Id eqId);
         
      protected:

      private:
         void linearStreamfunction(DecoupledZSparse& op, const MHDFloat k, const PhysicalNames::Id fieldId);
         void linearVelocityZ(DecoupledZSparse& op, const MHDFloat k, const PhysicalNames::Id fieldId);
         void linearTemperature(DecoupledZSparse& op, const MHDFloat k, const PhysicalNames::Id fieldId);
         void boundaryStreamfunction(DecoupledZSparse& op, const MHDFloat k, const PhysicalNames::Id fieldId);
         void boundaryVelocityZ(DecoupledZSparse& op, const MHDFloat k, const PhysicalNames::Id fieldId);
         void boundaryTemperature(DecoupledZSparse& op, const MHDFloat k, const PhysicalNames::Id fieldId);
   };

}
}

#endif // BETA3DQGSYSTEM_HPP
