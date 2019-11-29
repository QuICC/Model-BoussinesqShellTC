/** 
 * @file PolySetup.hpp
 * @brief Implementation of the polynomial transform setup class
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
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

namespace QuICC {

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
          * @param mult       List of second dimension multipliers
          * @param padSize    Padding size
          */
         PolySetup(const int size, const int howmany, const int specSize, const std::vector<ArrayI>& fast, const ArrayI& slow, const ArrayI& mult, const int padSize);

         /**
          * @brief Empty destructor
          */
         ~PolySetup();

         /**
          * @brief Get the fast indexes
          */
         const std::vector<ArrayI>&  fast() const;

         /**
          * @brief Get the slow indexes
          */
         const ArrayI&  slow() const;

         /**
          * @brief Get the multipliers
          */
         const ArrayI&  mult() const;

         /**
          * @brief Get padding size
          */
         int padSize() const;
         
      protected:

      private:
         /**
          * @brief Storage for the fast indexes
          */
         std::vector<ArrayI> mFast;

         /**
          * @brief Storage for the slow indexes
          */
         ArrayI mSlow;

         /**
          * @brief Storage for the multipliers
          */
         ArrayI mMult;

         /**
          * @brief Padding size
          */
         int mPadSize;
   };

   /// Typedef for an smart reference counting pointer for a PolySetup
   typedef SharedPtrMacro<PolySetup>   SharedPolySetup;

}
}

#endif // POLYSETUP_HPP
