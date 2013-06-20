/** \file PolySetup.hpp
 *  \brief Implementation of the polynomial transform setup class
 */

#ifndef POLYSETUP_HPP
#define POLYSETUP_HPP

// Configuration includes
//
#include "SmartPointers/SharedPtrMacro.h"

// System includes
//

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Resolutions/TransformSetup.hpp"

namespace GeoMHDiSCC {

namespace Transform {

   /**
    * @brief Implementation of the polynomial transform setup class
    */ 
   class PolySetup: public TransformSetup
   {
      public:
         /**
          * @brief Constructor
          *
          * @param size       Size of the transform
          * @param howmany    Number of similar transforms
          * @param specSize   Spectral output size (i.e without the padding)
          * @param fast       List of fast varying indexes
          * @param slow       List of slow varying indexes
          */
         PolySetup(const int size, const int howmany, const int specSize, const std::vector<ArrayI>& fast, const ArrayI& slow);

         /**
          * @brief Empty destructor
          */
         ~PolySetup();
         
      protected:

      private:
         /**
          * @brief Initialise the setup
          */
         void init();

         /**
          * @brief Storage for the fast indexes
          */
         std::vector<ArrayI> mFast;

         /**
          * @brief Storage for the slow indexes
          */
         ArrayI mSlow;
   };

   /// Typedef for an smart reference counting pointer for a PolySetup
   typedef SharedPtrMacro<PolySetup>   SharedPolySetup;

}
}

#endif // POLYSETUP_HPP
