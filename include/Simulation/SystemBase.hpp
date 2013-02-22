/** \file SystemBase.hpp
 *  \brief The lowest level system building block
 */

#ifndef SYSTEMBASE_HPP
#define SYSTEMBASE_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "Base/Resolutions/Resolution.hpp"
#include "Simulation/General/ExecutionTimer.hpp"

namespace EPMPhoenix {

   /**
    * \brief The lowest level system building block
    *
    * This is the building block for all the simulations and tools
    */
   class SystemBase
   {
      public:
         /**
          * @brief Constructor
          */
         SystemBase();

         /**
          * @brief Simple empty destructor
          */
         virtual ~SystemBase() {};

         /**
          * @brief Get the resolution information
          */
         SharedResolution spRes() const;

      protected:
         /**
          * @brief General execution timer
          */
         ExecutionTimer mExecTimer;

         /**
          * @brief Shared pointer to the resolution object
          */
         SharedResolution  mspRes;

         /**
          * @brief Initialise the base system
          */
         void initSystem();

         /**
          * @brief Cleanup unused memory from the base system
          */
         void cleanupSystem();

      private:
   };

   inline SharedResolution SystemBase::spRes() const
   {
      return this->mspRes;
   }

}

#endif // SYSTEMBASE_HPP
