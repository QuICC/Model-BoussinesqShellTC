/** \file FftSetup.hpp
 *  \brief Implementation of the FFT setup class
 */

#ifndef FFTSETUP_HPP
#define FFTSETUP_HPP

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

namespace GeoMHDiSCC {

namespace Transform {

   /**
    * @brief Implementation of the FFT setup class
    */ 
   class FftSetup
   {
      public:
         /**
          * @brief Constructor
          *
          * @param size       Size of the transform
          * @param howmany    Number of similar transforms
          * @param specSize   Spectral output size (i.e without the padding)
          * @param isMixed    FFT involves real and complex values (ie C2R and R2C)
          */
         FftSetup(const int size, const int howmany, const int specSize, const bool isMixed);

         /**
          * @brief Empty destructor
          */
         ~FftSetup() {};

         /**
          * @brief Does setup define a mixed real/complex transform?
          */
         bool isMixed() const;

         /**
          * @brief Set the transform scaling factor
          */
         void setScale(const MHDFloat scale);

         /**
          * @brief Get the size of the transform
          */
         int fwdSize() const;

         /**
          * @brief Get the size of the transform
          */
         int bwdSize() const;

         /**
          * @brief Get the number of similar transforms
          */
         int howmany() const;

         /**
          * @brief Get the spectral size of the transform
          */
         int specSize() const;

         /**
          * @brief Get the size of the padding
          */
         int padSize() const;

         /**
          * @brief Get scaling factor
          */
         MHDFloat scale() const;
         
      protected:

      private:
         /**
          * @brief Size of the forward transform
          */
         int mFwdSize;

         /**
          * @brief Size of the backward transform
          */
         int mBwdSize;

         /**
          * @brief Number of similar transforms
          */
         int mHowmany;

         /**
          * @brief Spectral size
          */
         int mSpecSize;

         /**
          * @brief Transform is of C2R/R2C type?
          */
         bool mIsMixed;

         /**
          * @brief Storage for the scale factor
          */
         MHDFloat mScale;
   };

   /// Typedef for an smart reference counting pointer for a FftSetup
   typedef SharedPtrMacro<FftSetup>   SharedFftSetup;

}
}

#endif // FFTSETUP_HPP
