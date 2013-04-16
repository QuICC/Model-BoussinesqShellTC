/** \file SerialDebugger.hpp
 *  \brief Implementation of a very simple serial debugger
 */

#ifndef SERIALDEBUGGER_HPP
#define SERIALDEBUGGER_HPP

// System includes
//
#include <string>

// External includes
//

// Project includes
//

namespace GeoMHDiSCC {

namespace Debug {

   /**
    * \brief Implementation of a very simple serial debugger
    */
   class SerialDebugger
   {
      public:
         /**
          * @brief Debug message when entering function
          *
          * @param msg Message to print
          * @param tabs Number of tab characters
          */
         static void enter(const std::string& msg, const int tabs);

         /**
          * @brief Debug message when leaving function
          *
          * @param msg Message to print
          * @param tabs Number of tab characters
          */
         static void leave(const std::string& msg, const int tabs);
         
      protected:

      private:
         /**
          * @brief Constructor
          */
         SerialDebugger();

         /**
          * @brief Destructor
          */
         ~SerialDebugger();
   };

}
}

#endif // SERIALDEBUGGER_HPP
