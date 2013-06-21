/** \file TransformSetup.hpp
 *  \brief Implementation of base class for a generalized transform setup
 */

#ifndef TRANSFORMSETUP_HPP
#define TRANSFORMSETUP_HPP

// Configuration includes
//
#include "SmartPointers/SharedPtrMacro.h"

// System includes
//

// External includes
//

// Project includes
//

namespace GeoMHDiSCC {

namespace Transform {

   /**
    * @brief Implementation of base class for a generalized transform setup
    */ 
   class TransformSetup
   {
      public:
         /**
          * @brief Constructor
          *
          * @param size       Size of the transform
          * @param howmany    Number of similar transforms
          * @param specSize   Spectral output size (i.e without the padding)
          */
         TransformSetup(const int size, const int howmany, const int specSize);

         /**
          * @brief Empty destructor
          */
         ~TransformSetup();

         /**
          * @brief Get the size of the transform
          */
         int fwdSize() const;

         /**
          * @brief Get the number of similar transforms
          */
         int howmany() const;

         /**
          * @brief Get the spectral size of the transform
          */
         int specSize() const;
         
      protected:
         /**
          * @brief Size of the forward transform
          */
         int mFwdSize;

         /**
          * @brief Number of similar transforms
          */
         int mHowmany;

         /**
          * @brief Spectral size
          */
         int mSpecSize;

      private:

   };

   /// Typedef for an smart reference counting pointer for a TransformSetup
   typedef SharedPtrMacro<TransformSetup>   SharedTransformSetup;

}
}

#endif // TRANSFORMSETUP_HPP
