/** \file Variable.hpp
 *  \brief Implementation of a generic simulation variable
 *
 *  \mhdBug Needs test
 */

#ifndef VARIABLE_HPP
#define VARIABLE_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "Resolutions/Resolution.hpp"
#include "Variables/VariableDomain.hpp"

namespace GeoMHDiSCC {

namespace Datatypes {

   /**
    * @brief This class implements the abstract concept of a simulation variable
    *
    * It is independent of it being a vector field or a scalar field.
    *
    * \tparam TVariable Type of the field
    */
   template <typename TVariable, int DOMAINS> class Variable : public VariableDomain<TVariable, DOMAINS>
   {
      public:
         /**
          * @brief Constructor
          *
          * @param spRes Resolution information
          */
         Variable(SharedResolution spRes);

         /**
          * @brief Simple empty destructor
          */
         virtual ~Variable();
         
      protected:

      private:
   };

   template <typename TVariable, int DOMAINS> Variable<TVariable,DOMAINS>::Variable(SharedResolution spRes)
      : VariableDomain<TVariable,DOMAINS>(spRes)
   {
   }

   template <typename TVariable, int DOMAINS> Variable<TVariable,DOMAINS>::~Variable()
   {
   }

}
}

#endif // VARIABLE_HPP
