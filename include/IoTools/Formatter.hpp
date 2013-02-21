/** \file Formatter.hpp
 *  \brief Implementation of a few useful output formating static functions
 */

#ifndef FORMATTER_HPP
#define FORMATTER_HPP

// Configuration includes
//

// System includes
//
#include <string>
#include <iostream>
#include <iomanip>

// External includes
//

// Project includes
//

namespace GeoMHDiSCC {

namespace IoTools {

   /**
    * \brief Implementation of a few useful output formating static functions
    */
   class Formatter
   {
      public:
         /**
          * @brief Base width for stdout output
          */
         static const int TEXT_WIDTH = 50;

         /**
          * @brief Output a newline
          */
         static void printNewline(std::ostream& stream);

         /**
          * @brief Output a line of characters
          */
         static void printLine(std::ostream& stream, char fill);

         /**
          * @brief Output a centered text with fill characters
          */
         static void printCentered(std::ostream& stream, const std::string& text, char fill = ' ', int base = -1);
         
      protected:

      private:
         /**
         * @brief Constructor 
         */
         Formatter();

         /**
         * @brief Destructor
         */
         ~Formatter() {};
   };

}
}

#endif // FORMATTER_HPP
