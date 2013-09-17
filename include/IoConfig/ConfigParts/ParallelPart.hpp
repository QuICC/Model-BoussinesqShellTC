/** 
 * @file ParallelPart.hpp 
 * @brief Implementation of the parallel part of the configuration file
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef PARALLELPART_HPP
#define PARALLELPART_HPP

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

namespace GeoMHDiSCC {

namespace IoConfig {

   /**
    * @brief Implementation of the parallel part of the configuration file
    */
   class ParallelPart: public IConfigurationPart
   {
      public:
         /**
          * @brief Constructor
          */
         ParallelPart();

         /**
          * @brief Destructor
          */
         virtual ~ParallelPart();

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

   /// Typedef for a shared pointer of a parallel part
   typedef SharedPtrMacro<ParallelPart> SharedParallelPart;

   /// Typedef for a const shared pointer of a parallel part
   typedef SharedPtrMacro<const ParallelPart> SharedCParallelPart;

}
}

#endif // PARALLELPART_HPP
