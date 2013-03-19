/** \file DiagnosticCoordinator.hpp
 *  \brief Coordinator for the diagnostics computations
 */

#ifndef DIAGNOSTICCOORDINATOR_HPP
#define DIAGNOSTICCOORDINATOR_HPP

// Configuration includes
//

// System includes
//
#include <string>

// External includes
//

// Project includes
//
#include "Enums/PhysicalNames.hpp"
#include "Diagnostics/IVelocityWrapper.hpp"
#include "TypeSelectors/VariableSelector.hpp"

namespace GeoMHDiSCC {

namespace Diagnostics {

   /**
    * @brief Coordinator for the diagnostics computations
    */
   class DiagnosticCoordinator
   {
      public:
         /**
          * @brief Constructor
          */
         DiagnosticCoordinator();

         /**
          * @brief Constructor
          */
         ~DiagnosticCoordinator();

         /**
          * @brief Initialise the coordinator
          *
          * @param mesh    Vector of grid values
          * @param scalars Map of shared scalar variables
          * @param vectors Map of shared vector variables
          */
         void init(const std::vector<Array>& mesh, const std::map<PhysicalNames::Id, Datatypes::SharedScalarVariableType>&  scalars, const std::map<PhysicalNames::Id, Datatypes::SharedVectorVariableType>&  vectors); 

         /**
          * @brief Compute the current CFL condition
          */
         void updateCfl();

         /**
          * @brief Compute the current kinetic energy
          */
         void updateKineticEnergy();

         /**
          * @brief Get CFL condition
          */
         MHDFloat cfl() const;

         /**
          * @brief Get kinetic energy condition
          */
         MHDFloat kineticEnergy() const;

      protected:

      private:
         /**
          * @brief Courant constant used for the CFL computation
          */
         const MHDFloat mcCourant;

         /**
          * @brief Current CFL condition
          */
         MHDFloat mCfl;

         /**
          * @brief Current kinetic energy
          */
         MHDFloat mKinetic;

         /**
          * @brief Minimal spacing between grid points
          */
         MHDFloat mMinSpacing;

         /**
          * @brief Shared pointer to a velocity field wrapper
          */
         SharedIVelocityWrapper  mspVelocityWrapper;
   };
}
}

#endif // DIAGNOSTICCOORDINATOR_HPP
