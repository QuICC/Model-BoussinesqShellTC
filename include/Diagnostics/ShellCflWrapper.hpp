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
#include "Diagnostics/IVectorWrapper.hpp"

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
          */
         ShellCflWrapper(const SharedIVectorWrapper spVelocity, const std::map<NonDimensional::Id,MHDFloat>& params);

         /**
          * @brief Constructor
          *
          * @param Velocity wrapper
          * @param Magnetic wrapper
          */
         ShellCflWrapper(const SharedIVectorWrapper spVelocity, const SharedIVectorWrapper spMagnetic, const std::map<NonDimensional::Id,MHDFloat>& params);

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
