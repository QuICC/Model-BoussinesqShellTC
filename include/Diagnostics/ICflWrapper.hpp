/**
 * @file ICflWrapper.hpp
 * @brief Interface for the CFL constraint
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef ICFLWRAPPER_HPP
#define ICFLWRAPPER_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "Diagnostics/ICflWrapper.hpp"
#include "Diagnostics/IVectorWrapper.hpp"

namespace QuICC {

namespace Diagnostics {

   /**
    * @brief Interface for the CFL constraint
    */
   class ICflWrapper
   {
      public:
         /**
          * @brief Constructor
          *
          * @param Velocity wrapper
          */
         ICflWrapper(const SharedIVectorWrapper spVelocity);

         /**
          * @brief Constructor
          *
          * @param Velocity wrapper
          * @param Magnetic wrapper
          */
         ICflWrapper(const SharedIVectorWrapper spVelocity, const SharedIVectorWrapper spMagnetic);

         /**
          * @brief Constructor
          */
         ~ICflWrapper();

         /**
          * @brief Initialize wrapper
          */
         virtual void init(const std::vector<Array>& mesh) = 0;

         /**
          * @brief Get initial CFL constraint
          */
         virtual MHDFloat initialCfl() const = 0;

         /**
          * @brief Get CFL constraint
          */
         virtual MHDFloat cfl() const = 0;

      protected:
         /**
          * @brief Shared velocity wrapper
          */
         SharedIVectorWrapper mspVelocity;

         /**
          * @brief Shared magnetic wrapper
          */
         SharedIVectorWrapper mspMagnetic;

      private:
   };

   /// Typedef for a shared ICflWrapper
   typedef SharedPtrMacro<ICflWrapper> SharedICflWrapper;
}
}

#endif // ICFLWRAPPER_HPP
