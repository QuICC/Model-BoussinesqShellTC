/**
 * @file Exception.hpp
 * @brief Simple exception class 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 * @version 0.9.0
 * @date 2013-09-11
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
         explicit Exception(const std::string& msg);

         /**
         * @brief Simple empty destructor
         */
         virtual ~Exception() throw();
         
      protected:

      private:
   };

}

#endif // EXCEPTION_HPP
