/** \file Exception.hpp
 *  \brief Simple exception class.
 */

#ifndef EXCEPTION_HPP
#define EXCEPTION_HPP

// System includes
//
#include <stdexcept>

// Project includes
//

namespace GeoMHDiSCC {

   /**
    * @brief Simple exception class
    *
    * This is a very simple exception class. It allows to give a location and a message string.
    * These can then be recovered in a catch block.
    */
   class Exception: public std::runtime_error
   {
      public:
         /**
         * @brief Constructs the exception based on the given strings.
         *
         * \param msg Error message
         */
         Exception(const std::string& msg);

         /**
         * @brief Simple empty destructor
         */
         virtual ~Exception() {};
         
      protected:

      private:
   };

}

#endif // EXCEPTION_HPP
