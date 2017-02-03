/** 
 * @file RunPart.hpp 
 * @brief Implementation of the run part of the configuration file
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef RUNPART_HPP
#define RUNPART_HPP

// Configuration includes
//

// System includes
//
#include <string>

// External includes
//

// Project includes
//
#include "IoConfig/ConfigParts/IConfigurationPart.hpp"

namespace QuICC {

namespace IoConfig {

   /**
    * @brief Implementation of the run part of the configuration file
    */
   class RunPart: public IConfigurationPart
   {
      public:
         /**
          * @brief Constructor
          */
         RunPart();

         /**
          * @brief Destructor
          */
         virtual ~RunPart();

         /**
          * @brief Check compatibility of data
          *
          * \mhdBug Check is not yet implemented
          */
         virtual void checkData();
         
      protected:
         /**
          * @brief Tag name of the parent node
          */
         static const std::string  PARENTTAG;

         /**
          * @brief Initialise component
          */
         void init();

      private:
   };

   /// Typedef for a shared pointer of a run part
   typedef SharedPtrMacro<RunPart> SharedRunPart;

   /// Typedef for a const shared pointer of a run part
   typedef SharedPtrMacro<const RunPart> SharedCRunPart;

}
}

#endif // RUNPART_HPP
