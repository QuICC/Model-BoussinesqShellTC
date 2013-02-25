/** \file CoreResolution.hpp
 *  \brief Definition of a resolution object for a single CPU
 */

#ifndef CORERESOLUTION_HPP
#define CORERESOLUTION_HPP

// Configuration includes
//
#include "SmartPointers/SharedPtrMacro.h"

// System includes
//
#include <vector>

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Resolutions/TransformResolution.hpp"

namespace GeoMHDiSCC {

   /**
    * @brief Definition of a resolution object for a single CPU
    */
   class CoreResolution
   {
      public:
         /**
          * @brief Constructor
          *
          * @param transformRes Resolution object for the transforms
          */
         explicit CoreResolution(const std::vector<SharedTransformResolution>& transformRes);

         /**
          * @brief Empty Destructor
          */
         ~CoreResolution();

         /**
          * @brief Get transform resolution for corresponding dimension
          *
          * @param i Index of the dimension
          */
         SharedCTransformResolution dim(const int i) const;

         /**
          * @brief Get number of dimensions
          */
         int nDim() const;

      protected:

      private:
         /**
          * @brief Resolution information for all transforms
          */
         std::vector<SharedTransformResolution>   mTransforms;
   };

   /// Typedef for a shared pointer to a CoreResolution object
   typedef SharedPtrMacro<CoreResolution>   SharedCoreResolution;

   /// Typedef for a shared pointer to a const CoreResolution object
   typedef SharedPtrMacro<const CoreResolution>   SharedCCoreResolution;

}

#endif // CORERESOLUTION_HPP
