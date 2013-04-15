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
          * @param masg Message to print
          */
         static void enter(const std::string& msg);

         /**
          * @brief Debug message when leaving function
          *
          * @param masg Message to print
          */
         static void leave(const std::string& msg);
         
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
