/**
 * @file SphereCflWrapper.hpp
 * @brief CFL constraint in a full sphere 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef QUICC_DIAGNOSTICS_SPHERECFLWRAPPER_HPP
#define QUICC_DIAGNOSTICS_SPHERECFLWRAPPER_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "Diagnostics/ISphericalCflWrapper.hpp"
#include "Diagnostics/IVelocityWrapper.hpp"

namespace QuICC {

namespace Diagnostics {

   /**
    * @brief CFL constraint in a full sphere
    */
   class SphereCflWrapper: public ISphericalCflWrapper
   {
      public:
         /**
          * @brief Constructor
          *
          * @param Velocity wrapper
          * @param Vector of physical space grid
          */
         SphereCflWrapper(const SharedIVelocityWrapper spVelocity);

         /**
          * @brief Constructor
          */
         ~SphereCflWrapper();

      protected:

      private:
         /**
          * @brief Get effective max harmonic degree L
          */
         virtual MHDFloat effectiveMaxL(const MHDFloat r) const;

         /**
          * @brief Compute approximation to first Jacobi root for CFL
          */
         MHDFloat jacobiRoot(const MHDFloat l) const;
   };

   /// Typedef for a shared SphereCflWrapper
   typedef SharedPtrMacro<SphereCflWrapper> SharedSphereCflWrapper;
}
}

#endif // QUICC_DIAGNOSTICS_SPHERECFLWRAPPER_HPP
