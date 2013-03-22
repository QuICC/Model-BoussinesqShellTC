/** \file SerialFramework.hpp
 *  \brief Implementation of the serial framework
 */

#ifndef SERIALFRAMEWORK_HPP
#define SERIALFRAMEWORK_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Framework/FrameworkBase.hpp"

namespace GeoMHDiSCC {

   /**
    * @brief This class defines the framework for a serial code
    */
   class SerialFramework: public FrameworkBase
   {
      public:
         /**
          * @brief Initialise the Serial system
          */
         static void init();

         /**
          * @brief Setup the Serial system
          */
         static void setup(const int nodes);

         /**
          * @brief Synchronise
          */
         static void synchronize();
  
         /**
          * @brief Finalise the Serial system
          */
         static void finalize();
         
      protected:

      private:
         /**
          * @brief Constructor
          */
         SerialFramework();

         /**
          * @brief Destructor
          */
         ~SerialFramework();
   };

}

#endif // SERIALFRAMEWORK_HPP
