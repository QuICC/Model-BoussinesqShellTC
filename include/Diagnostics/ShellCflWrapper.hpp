/**
 * @file ShellCflWrapper.hpp
 * @brief CFL constraint in a spherical shell 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef QUICC_DIAGNOSTICS_SHELLCFLWRAPPER_HPP
#define QUICC_DIAGNOSTICS_SHELLCFLWRAPPER_HPP

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
    * @brief CFL constraint in a spherical Shell
    */
   class ShellCflWrapper: public ISphericalCflWrapper
   {
      public:
         /**
          * @brief Constructor
          *
          * @param Velocity wrapper
          * @param Vector of physical space grid
          */
         ShellCflWrapper(const SharedIVelocityWrapper spVelocity);

         /**
          * @brief Constructor
          */
         ~ShellCflWrapper();

      protected:

      private:
         /**
          * @brief Get effective max harmonic degree L
          */
         virtual MHDFloat effectiveMaxL(const MHDFloat r) const;
   };

   /// Typedef for a shared ShellCflWrapper
   typedef SharedPtrMacro<ShellCflWrapper> SharedShellCflWrapper;
}
}

#endif // QUICC_DIAGNOSTICS_SHELLCFLWRAPPER_HPP
