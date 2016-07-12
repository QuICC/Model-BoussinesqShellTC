/**
 * @file FftSetup.hpp
 * @brief Implementation of the FFT setup class 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
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
#include "Resolutions/TransformSetup.hpp"

namespace GeoMHDiSCC {

namespace Transform {

   /**
    * @brief Implementation of the FFT setup class
    */ 
   class FftSetup: public TransformSetup
   {
      public:
         /**
          * List of possible size combinations
          */
         enum Type {
            /// Real <-> real
            REAL,
            /// Complex <-> Complex
            COMPLEX,
            /// Complex -> real, real -> complex sizes
            MIXED,
            /// Complex -> complex but componentwise
            COMPONENT
         };

         /**
          * @brief Constructor
          *
          * @param size       Size of the transform
          * @param howmany    Number of similar transforms
          * @param specSize   Spectral output size (i.e without the padding)
          * @param type       Type of the transform involved
          */
         FftSetup(const int size, const int howmany, const int specSize, const FftSetup::Type type);

         /**
          * @brief Constructor
          *
          * @param size       Size of the transform
          * @param howmany    Number of similar transforms slip in even and odd
          * @param evenBlocks Distribution of the even blocks in full data
          * @param oddBlocks  Distribution of the odd blocks in full data
          * @param specSize   Spectral output size (i.e without the padding)
          * @param type       Type of the transform involved
          */
         FftSetup(const int size, const ArrayI& howmany, const MatrixI& evenBlocks, const MatrixI& oddBlocks, const int specSize, const FftSetup::Type type);

         /**
          * @brief Constructor
          *
          * @param size          Size of the transform
          * @param howmany       Number of similar transforms slip in even and odd
          * @param idBlocks      Block sizes for 3D index
          * @param specSize      Spectral output size (i.e without the padding)
          * @param type          Type of the transform involved
          */
         FftSetup(const int size, const int howmany, const MatrixI& idBlocks, const int specSize, const FftSetup::Type type);

         /**
          * @brief Empty destructor
          */
         ~FftSetup();

         /**
          * @brief Does setup define a mixed real/complex transform?
          */
         FftSetup::Type type() const;

         /**
          * @brief Set the transform scaling factor
          */
         void setScale(const MHDFloat scale);

         /**
          * @brief Set the box size scaling factor
          */
         void setBoxScale(const MHDFloat boxScale);

         /**
          * @brief Get the size of the transform
          */
         int bwdSize() const;

         /**
          * @brief Get the size of the padding
          */
         int padSize() const;

         /**
          * @brief Get scaling factor
          */
         MHDFloat scale() const;

         /**
          * @brief Get box size scaling factor
          */
         MHDFloat boxScale() const;

         /**
          * @brief Get the even block distribution
          */
         const MatrixI& evenBlocks() const; 

         /**
          * @brief Get the odd block distribution
          */
         const MatrixI& oddBlocks() const; 

         /**
          * @brief Get the block distribution with 3D index
          */
         const MatrixI& idBlocks() const; 
         
      protected:

      private:
         /**
          * @brief Set backward size
          */
         void setBackwardSize();

         /**
          * @brief Size of the backward transform
          */
         int mBwdSize;

         /**
          * @brief Transform is of C2R/R2C type?
          */
         FftSetup::Type mType;

         /**
          * @brief Storage for the scale factor
          */
         MHDFloat mScale;

         /**
          * @brief Storage for the box scale factor for derivatives
          */
         MHDFloat mBoxScale;

         /**
          * @brief Storage for the even blocks distribution description
          */
         MatrixI mEvenBlocks;

         /**
          * @brief Storage for the even blocks distribution description
          */
         MatrixI mOddBlocks;

         /**
          * @brief Storage for the block sizes with 3D index
          */
         MatrixI mIdBlocks;
   };

   /// Typedef for an smart reference counting pointer for a FftSetup
   typedef SharedPtrMacro<FftSetup>   SharedFftSetup;

}
}

#endif // FFTSETUP_HPP
