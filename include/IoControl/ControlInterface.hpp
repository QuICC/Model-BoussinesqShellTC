/** \file ControlInterface.hpp
 *  \brief Implementation of a external runtime control file
 */

#ifndef CONTROLINTERACE_HPP
#define CONTROLINTERACE_HPP

// System includes
//
#include <fstream>

// External includes
//

// Project includes
//
#include "Enums/RuntimeStatus.hpp"
#include "IoControl/ControlFile.hpp"

namespace GeoMHDiSCC {

namespace IoControl {

   /**
    * @brief Implementation of an external runtime control file
    */
   class ControlInterface: public ControlFile
   {
      public:
         /**
         * @brief Constructor
         */
         ControlInterface();

         /**
         * @brief Destructor
         */
         ~ControlInterface();

         /**
          * @brief Read input or check existance
          */
         void read();

         /**
          * @brief Current runtime status?
          */
         RuntimeStatus::Id status() const;
         
      protected:

      private:
         /**
          * @brief Handle to the file
          */
         std::fstream mFile;

         /**
          * @brief Current runtime status
          */
         RuntimeStatus::Id mStatus;

         /**
          * @brief Initialise control file
          */
         void init();

         /**
          * @brief Finalise the file
          */
         void finalize();

         /**
          * @brief Create the file
          */
         void create();

         /**
          * @brief Open the file
          */
         void open();

         /**
          * @brief Close the file
          */
         void close();

         /**
          * @brief Delete the file
          */
         void deleteFile();

         /**
          * @brief Update status
          */
         void update();
   };
}
}

#endif // CONTROLINTERACE_HPP
