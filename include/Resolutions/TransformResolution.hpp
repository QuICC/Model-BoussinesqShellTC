/** 
 * @file TransformResolution.hpp
 * @brief Definition of the required resolution information for a transform
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef TRANSFORMRESOLUTION_HPP
#define TRANSFORMRESOLUTION_HPP

// Configuration includes
//
#include "StaticAsserts/StaticAssert.hpp"

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
#include "TypeSelectors/ScalarSelector.hpp"

namespace QuICC {

   /**
    * @brief Definition of the required resolution information for a transform
    *
    * \mhdBug The resolution with dimension only is not safe and will segfault if idx methods are called
    */
   class TransformResolution
   {
      public:
         /**
          * @brief Constructor for full resolution
          *
          * @param fwd Set of indexes for 1D (forward direction)
          * @param bwd Set of indexes for 1D (backward direction)
          * @param idx2D Set of indexes for 2D
          * @param idx3D Set of indexes for 3D
          */
         TransformResolution(const std::vector<ArrayI>& fwd, const std::vector<ArrayI>& bwd, const std::vector<ArrayI>& idx2D, const ArrayI& idx3D);

         /**
          * @brief Empty Destructor
          */
         ~TransformResolution();

         /**
          * @brief Get dimension
          *
          * @param k Index of the last dimension
          */
         template <Dimensions::Data::Id TDim> int dim(const int k) const;

         /**
          * @brief Get dimension
          */
         template <Dimensions::Data::Id TDim> int dim() const;

         /**
          * @brief Get index mapping to full resolution
          *
          * @param i Index of dimension
          * @param k Index of the last dimension
          */
         template <Dimensions::Data::Id TDim> int idx(const int i, const int k) const;

         /**
          * @brief Get index mapping to full resolution
          *
          * @param i Index of dimension
          */
         template <Dimensions::Data::Id TDim> int idx(const int i) const;

         /**
          * @brief Get mode index mapping to full resolution indexes
          *
          * @param i Mode index
          */
         ArrayI mode(const int i) const;

         /**
          * @brief Clear the indexes
          */
         void clearIndexes();

      protected:

      private:
         /**
          * @brief Initialise the dimensions from the indexes
          */
         void initDimensions();

         /**
          * @brief Set of indexes describing the first dimensions (forward direction)
          */
         std::vector<ArrayI>   mFwd;

         /**
          * @brief Set of indexes describing the first dimensions (backward direction)
          */
         std::vector<ArrayI>   mBwd;

         /**
          * @brief Set of indexes describing the second dimensions
          */
         std::vector<ArrayI>   mIdx2D;

         /**
          * @brief Set of indexes describing the third dimensions
          */
         ArrayI   mIdx3D;

         /**
          * @brief Dimension of the first dimension (forward direction)
          */
         ArrayI   mDimF1D;

         /**
          * @brief Dimension of the first dimension (backward direction)
          */
         ArrayI   mDimB1D;

         /**
          * @brief Dimension of the second dimension
          */
         ArrayI   mDim2D;

         /**
          * @brief Dimension of the third dimension
          */
         int mDim3D;
   };

   template <Dimensions::Data::Id TID> int TransformResolution::dim(const int k) const
   {
      // This should never by used
      Debug::StaticAssert< TID == -99 >();

      return -1;
   }

   template <Dimensions::Data::Id TID> int TransformResolution::dim() const
   {
      // This should never by used
      Debug::StaticAssert< TID == -99 >();

      return -1;
   }

   template <Dimensions::Data::Id TID> int TransformResolution::idx(const int i, const int k) const
   {
      // This should never by used
      Debug::StaticAssert< TID == -99 >();

      return -1;
   }

   template <Dimensions::Data::Id TID> int TransformResolution::idx(const int i) const
   {
      // This should never by used
      Debug::StaticAssert< TID == -99 >();

      return -1;
   }

   /// Specialised for DATF1D
   template  <> int TransformResolution::dim<Dimensions::Data::DATF1D>(const int k) const;
   /// Specialised for DATB1D
   template  <> int TransformResolution::dim<Dimensions::Data::DATB1D>(const int k) const;
   /// Specialised for DAT2D
   template  <> int TransformResolution::dim<Dimensions::Data::DAT2D>(const int k) const;
   /// Specialised for DATF1D
   template  <> int TransformResolution::dim<Dimensions::Data::DATF1D>() const;
   /// Specialised for DATB1D
   template  <> int TransformResolution::dim<Dimensions::Data::DATB1D>() const;
   /// Specialised for DAT2D
   template  <> int TransformResolution::dim<Dimensions::Data::DAT2D>() const;
   /// Specialised for DAT3D
   template  <> int TransformResolution::dim<Dimensions::Data::DAT3D>() const;

   /// Specialised for DATF1D
   template  <> int TransformResolution::idx<Dimensions::Data::DATF1D>(const int i, const int k) const;
   /// Specialised for DATB1D
   template  <> int TransformResolution::idx<Dimensions::Data::DATB1D>(const int i, const int k) const;
   /// Specialised for DAT2D
   template  <> int TransformResolution::idx<Dimensions::Data::DAT2D>(const int i, const int k) const;
   /// Specialised for DATF1D
   template  <> int TransformResolution::idx<Dimensions::Data::DATF1D>(const int i) const;
   /// Specialised for DATB1D
   template  <> int TransformResolution::idx<Dimensions::Data::DATB1D>(const int i) const;
   /// Specialised for DAT2D
   template  <> int TransformResolution::idx<Dimensions::Data::DAT2D>(const int i) const;
   /// Specialised for DAT3D
   template  <> int TransformResolution::idx<Dimensions::Data::DAT3D>(const int i) const;

   /// Typedef for a shared pointer to a TransformResolution object
   typedef SharedPtrMacro<TransformResolution>   SharedTransformResolution;

   /// Typedef for a shared pointer to a const TransformResolution object
   typedef SharedPtrMacro<const TransformResolution>   SharedCTransformResolution;

}

#endif // TRANSFORMRESOLUTION_HPP
