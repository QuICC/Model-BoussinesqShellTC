/**
 * @file Arithmetics.hpp
 * @brief Definition of some useful enums for basic arithmetic operations 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
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

namespace QuICC {

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
         /// Set new negative value
         SETNEG,
         /// Add to value
         ADD,
         /// Substract from value
         SUB,
         /// Do nothing
         NONE,
      };
   };
}

#endif // ARITHMETICS_HPP
