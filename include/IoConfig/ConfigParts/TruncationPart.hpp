/** \file TruncationPart.hpp 
 *  \brief Implementation of the truncation part of the configuration file
 *
 *  \mhdBug Needs test
 */

#ifndef TRUNCATIONPART_HPP
#define TRUNCATIONPART_HPP

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
    * @brief Implementation of the truncation part of the configuration file
    *
    * \mhdTodo Strides have been desactivated (need first a clear idea how to use them)
    */
   class TruncationPart: public IConfigurationPart
   {
      public:
         /**
          * @brief Constructor
          *
          * @param dim Dimensionality of truncation
          */
         explicit TruncationPart(const int dim, const std::vector<bool>& isPeriodicBox);

         /**
          * @brief Destructor
          */
         virtual ~TruncationPart();

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
          *
          * @param dim Dimensionality of truncation
          */
         void init(const int dim, const std::vector<bool>& isPeriodicBox);

      private:

   };

   /// Typedef for a shared pointer of a truncation part
   typedef SharedPtrMacro<TruncationPart> SharedTruncationPart;

   /// Typedef for a const shared pointer of a truncation part
   typedef SharedPtrMacro<const TruncationPart> SharedCTruncationPart;

}
}

#endif // TRUNCATIONPART_HPP
