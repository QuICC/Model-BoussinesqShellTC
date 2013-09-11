/** 
 * @file VariableBase.hpp
 * @brief Base of the implementation of the variables
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 * @version 0.9.0
 * @date 2013-09-11
 */

#ifndef VARIABLEBASE_HPP
#define VARIABLEBASE_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "Resolutions/Resolution.hpp"

namespace GeoMHDiSCC {

namespace Datatypes {

   /**
    * @brief Base of the implementation of the variables
    */
   class VariableBase
   {
      public:
         /**
         * @brief Construct the general shared information for a physical variable
         *
         * @param spRes Resolution information
         */
         VariableBase(SharedResolution spRes);

         /**
         * @brief Destructor
         */
         virtual ~VariableBase();

         /**
          * @brief Get resolution information
          */
         const SharedResolution  spRes() const;

     #ifdef GEOMHDISCC_STORAGEPROFILE
         /**
         * @brief Get the memory requirements
         */
         MHDFloat requiredStorage() const;
     #endif // GEOMHDISCC_STORAGEPROFILE
         
      protected:

      private:
         /**
          * @brief Pointer to resolution information
          */
         SharedResolution   mspRes;
   };

   inline const SharedResolution VariableBase::spRes() const
   {
      return this->mspRes;
   }

}
}

#endif // VARIABLEBASE_HPP
