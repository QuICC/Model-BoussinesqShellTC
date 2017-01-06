/**
 * @file EquationParameters.hpp
 * @brief Definition of the non-dimensional parameters
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef EQUATIONPARAMETERS_HPP
#define EQUATIONPARAMETERS_HPP

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

namespace QuICC {

namespace Equations {

   /**
    * @brief This class provides the constant coefficients apparearing in the equations
    */
   class EquationParameters
   {
      public:
         /**
          * @brief Constructor
          */
         EquationParameters();

         /**
          * @brief Destructor
          */
         virtual ~EquationParameters();

         /**
          * @brief Get the IDs of the equation parameters
          */
         std::vector<NonDimensional::Id>  ids() const;

         /**
          * @brief Get the names of the equation parameters
          */
         std::vector<std::string>  names() const;

         /**
          * @brief Initialise the values from given parameters
          *
          * @param parameters Parameter values read from configuration
          */
         void init(const std::map<std::string, MHDFloat>& parameters);

         /**
          * @brief Get a nondimensional parameter
          *
          * @param name Name of the parameter
          */
         MHDFloat nd(NonDimensional::Id name) const;

         /**
          * @brief Get full map of parameters
          */
         const std::map<NonDimensional::Id,MHDFloat>& map() const;

         /**
          * @brief Update full map of parameters
          */
         std::map<NonDimensional::Id,MHDFloat>& rMap();
         
      protected:
         /**
          * @brief Storage for the nondimensional parameters
          */
         std::map<NonDimensional::Id, MHDFloat>   mND;

      private:
   };

   /// Typedef for a shared pointer to an EquationParameters object
   typedef SharedPtrMacro<EquationParameters>   SharedEquationParameters;
}
}

#endif // EQUATIONPARAMETERS_HPP
