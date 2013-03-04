/** \file PrRaXGParameters.hpp
 *  \brief Definition of the \f$Pr\f$, \f$Ra\f$, \f$\chi\f$, \f$\Gamma\f$ based parameters
 *
 *  \mhdBug Needs test
 */

#ifndef PRRAXGPARAMETERS_HPP
#define PRRAXGPARAMETERS_HPP

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
#include "Equations/IEquationParameters.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   /**
    * @brief This class provides the constant coefficients apparearing in the equations
    */
   class PrRaXGParameters: public IEquationParameters
   {
      public:
         /**
          * @brief Constructor
          */
         PrRaXGParameters();

         /**
          * @brief Destructor
          */
         virtual ~PrRaXGParameters();

         /**
          * @brief Get the names of the equation parameters
          */
         virtual std::vector<std::string>  names() const;

         /**
          * @brief Initialise the values from given parameters
          *
          * @param parameters Parameter values read from configuration
          */
         virtual void init(const std::map<std::string, MHDFloat>& parameters);
         
      protected:

      private:
   };

   /// Typedef for a shared pointer to an PrRaXGParameters object
   typedef SharedPtrMacro<PrRaXGParameters>   SharedPrRaXGParameters;
}
}

#endif // PRRAXGPARAMETERS_HPP
