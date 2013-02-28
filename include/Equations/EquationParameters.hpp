/** \file EquationParameters.hpp
 *  \brief Definition of the non-dimensional parameters
 *
 *  \mhdBug Needs test
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
//#include "Simulation/Enums/ParameterNames.hpp"

namespace GeoMHDiSCC {

   /**
    * @brief This class provides the constant coefficients apparearing in the equations
    *
    * \mhdBug This implementation is very buggy
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
         virtual ~EquationParameters() {};
//
//         /**
//          * @brief Get the names of the equation parameters
//          */
//         virtual std::vector<std::string>  names() const = 0;
//
//         /**
//          * @brief Initialise the values from given parameters
//          *
//          * @param parameters Parameter values read from configuration
//          */
//         virtual void init(const std::map<std::string, EPMFloat>& parameters) = 0;
//
//         /**
//          * @brief Get a nondimensional parameter
//          *
//          * @param name Name of the parameter
//          */
//         EPMFloat nd(ParameterNames::NonDimensional::Id name) const;
//
//         /**
//          * @brief Get an equation factor
//          *
//          * @param name Name of the factor
//          */
//         EPMFloat eq(ParameterNames::Equation::Id name) const;
//
//         /**
//          * @brief Get an energy factor
//          *
//          * @param name Name of the factor
//          */
//         EPMFloat e(ParameterNames::Energy::Id name) const;
//
//         /**
//          * @brief Get a global CFL condition
//          *
//          * @param name Name of the condition
//          */
//         EPMFloat cfl(ParameterNames::Cfl::Id name) const;
//
//         /**
//          * @brief Get a field scaling factor
//          *
//          * @param name Name of the scaling factor
//          */
//         EPMFloat s(ParameterNames::Scaling::Id name) const;
//
//         /**
//          * @brief Get an extra factor
//          *
//          * @param name Name of the factor
//          */
//         EPMFloat x(ParameterNames::Extra::Id name) const;
         
      protected:
//         /**
//          * @brief Storage for the nondimensional parameters
//          */
//         std::map<ParameterNames::NonDimensional::Id, EPMFloat>   mND;
//
//         /**
//          * @brief Storage for the equation factors
//          */
//         std::map<ParameterNames::Equation::Id, EPMFloat>   mEq;
//
//         /**
//          * @brief Storage for the energy factors
//          */
//         std::map<ParameterNames::Energy::Id, EPMFloat>   mE;
//
//         /**
//          * @brief Storage for the CFL conditions
//          */
//         std::map<ParameterNames::Cfl::Id, EPMFloat>   mCfl;
//
//         /**
//          * @brief Storage for the scaling factors
//          */
//         std::map<ParameterNames::Scaling::Id, EPMFloat>   mS;
//
//         /**
//          * @brief Storage for the extra factors
//          */
//         std::map<ParameterNames::Extra::Id, EPMFloat>   mX;

      private:
   };

//   inline EPMFloat EquationParameters::nd(ParameterNames::NonDimensional::Id name) const
//   {
//      return this->mND.at(name);
//   }
//
//   inline EPMFloat EquationParameters::eq(ParameterNames::Equation::Id name) const
//   {
//      return this->mEq.at(name);
//   }
//
//   inline EPMFloat EquationParameters::e(ParameterNames::Energy::Id name) const
//   {
//      return this->mE.at(name);
//   }
//
//   inline EPMFloat EquationParameters::cfl(ParameterNames::Cfl::Id name) const
//   {
//      return this->mCfl.at(name);
//   }
//
//   inline EPMFloat EquationParameters::s(ParameterNames::Scaling::Id name) const
//   {
//      return this->mS.at(name);
//   }
//
//   inline EPMFloat EquationParameters::x(ParameterNames::Extra::Id name) const
//   {
//      return this->mX.at(name);
//   }

   /// Typedef for a shared pointer to an EquationParameters object
   typedef SharedPtrMacro<EquationParameters>   SharedEquationParameters;
}

#endif // EQUATIONPARAMETERS_HPP
