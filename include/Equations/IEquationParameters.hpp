/** \file IEquationParameters.hpp
 *  \brief Definition of the non-dimensional parameters
 */

#ifndef IEQUATIONPARAMETERS_HPP
#define IEQUATIONPARAMETERS_HPP

// Configuration includes
//
#include "SmartPointers/SharedPtrMacro.h"

// System includes
//
#include <map>

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Enums/NonDimensional.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   /**
    * @brief This class provides the constant coefficients apparearing in the equations
    */
   class IEquationParameters
   {
      public:
         /**
          * @brief Constructor
          */
         IEquationParameters();

         /**
          * @brief Destructor
          */
         virtual ~IEquationParameters();

         /**
          * @brief Initialise the values from given parameters
          *
          * @param parameters Parameter values read from configuration
          */
         virtual void init(const std::map<std::string, MHDFloat>& parameters) = 0;

         /**
          * @brief Get a nondimensional parameter
          *
          * @param name Name of the parameter
          */
         MHDFloat nd(NonDimensional::Id name) const;
         
      protected:
         /**
          * @brief Storage for the nondimensional parameters
          */
         std::map<NonDimensional::Id, MHDFloat>   mND;

      private:
   };

   /// Typedef for a shared pointer to an IEquationParameters object
   typedef SharedPtrMacro<IEquationParameters>   SharedIEquationParameters;
}
}

#endif // IEQUATIONPARAMETERS_HPP