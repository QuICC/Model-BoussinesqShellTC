/**
 * @file Arithmetics.hpp
 * @brief Definition of some useful enums for basic arithmetic operations 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 * @version 0.9.0
 * @date 2013-09-11
 */

#ifndef ARITHMETICS_HPP
#define ARITHMETICS_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//

namespace GeoMHDiSCC {

   /**
    * @brief Simple struct to hold the arithmetic operations
    */
   struct Arithmetics {

      /**
      * @name Enum for basic arithmetics to use a template argument
      */
      enum Id {
         /// Set new value
         SET,
         /// Add to value
         ADD,
         /// Substract from value
         SUB
      };
   };
}

#endif // ARITHMETICS_HPP
