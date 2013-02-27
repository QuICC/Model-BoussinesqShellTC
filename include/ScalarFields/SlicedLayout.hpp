/** \file SlicedLayout.hpp
 *  \brief Base for a  scalar field with sliced data layout
 */

#ifndef SLICEDLAYOUT_HPP
#define SLICEDLAYOUT_HPP

// Debug includes
//
#include "StorageProfiler/StorageProfilerMacro.h"
#include "Exceptions/Exception.hpp"
#include "StaticAsserts/StaticAssert.hpp"

// Configuration includes
//

// System includes
//
#include <vector>

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Enums/Dimensions.hpp"

namespace GeoMHDiSCC {

namespace Datatypes {

   /**
    * @brief Base for a  scalar field with sliced data layout
    */
   template <typename TData, Dimensions::Type DIMENSION> class SlicedLayout
   {
      public:
         /// Typedef for the coefficient type
         typedef TData CoefficientType;

         /// Typedef for the storage type
         typedef Eigen::Matrix<TData, Eigen::Dynamic, Eigen::Dynamic> StorageType;

         /**
          * @brief Constructor
          */
         SlicedLayout();

         /**
          * @brief Destructor
          */
         virtual ~SlicedLayout();

         /**
          * @brief Get sliced field data
          */
         const std::vector<StorageType>& sliced() const;

         /**
          * @brief Set sliced field data
          */
         std::vector<StorageType>& rSliced();

         /**
          * @brief Get internal storage field data
          */
         const std::vector<StorageType>& data() const;

         /**
          * @brief Set internal storage field data
          */
         std::vector<StorageType>& rData();

         /**
          * @brief Get a 2D slice of the field
          *
          * @param k Index of the slice
          */
         const StorageType& slice(const int k) const;

         /**
          * @brief Set a 2D slice of the field
          *
          * @param k Index of the slice
          */
         StorageType& rSlice(const int k);

         /**
          * @brief Get a 1D profiled of the field
          *
          * @param j Index of the profile
          * @param k Index of the slice
          */
         typename StorageType::ConstColXpr profile(const int j, const int k = 0) const;

         /**
          * @brief Set a 1D profiled of the field
          *
          * @param j Index of the profile
          * @param k Index of the slice
          */
         typename StorageType::ColXpr rProfile(const int j, const int k = 0);

         /**
          * @brief Get a point value of the field
          *
          * @param i Index of the point
          * @param j Index of the profile
          * @param k Index of the slice
          */
         const TData& point(const int i, const int j = 0, const int k = 0) const;

         /**
          * @brief Set a point value of the field
          *
          * @param i Index of the point
          * @param j Index of the profile
          * @param k Index of the slice
          */
         TData& rPoint(const int i, const int j = 0, const int k = 0);

         /**
          * @brief Initialise field to zero
          */
         void initialiseZeros();

         /**
          * @brief Rescale the coefficients by a constant
          *
          * @param scale Scaling factor
          */
         void rescale(const MHDFloat scale);

     #ifdef GEOMHDISCC_STORAGEPROFILE
         /**
         * @brief Get the memory requirements
         */
         MHDFloat requiredStorage() const;
     #endif // GEOMHDISCC_STORAGEPROFILE
         
      protected:
         /**
          * @brief Field values
          */
         std::vector<StorageType>  mField;

         /**
          * @brief Initialise the storage
          *
          * @param dim1D Size in the first dimension
          * @param dim2D Size in the second dimension
          * @param dim3D Size in the third dimension
          */
         void initStorage(const ArrayI& dim1D, const ArrayI& dim2D, const int dim3D);

      private:
   };

   template <typename TData, Dimensions::Type DIMENSION> inline const std::vector<typename SlicedLayout<TData, DIMENSION>::StorageType>& SlicedLayout<TData, DIMENSION>::sliced() const
   {
      return this->mField;
   }

   template <typename TData, Dimensions::Type DIMENSION> inline std::vector<typename SlicedLayout<TData, DIMENSION>::StorageType>& SlicedLayout<TData, DIMENSION>::rSliced()
   {
      return this->mField;
   }

   template <typename TData, Dimensions::Type DIMENSION> inline const std::vector<typename SlicedLayout<TData, DIMENSION>::StorageType>& SlicedLayout<TData, DIMENSION>::data() const
   {
      return this->mField;
   }

   template <typename TData, Dimensions::Type DIMENSION> inline std::vector<typename SlicedLayout<TData, DIMENSION>::StorageType>& SlicedLayout<TData, DIMENSION>::rData()
   {
      return this->mField;
   }

   template <typename TData, Dimensions::Type DIMENSION> inline const typename SlicedLayout<TData, DIMENSION>::StorageType& SlicedLayout<TData, DIMENSION>::slice(const int k) const
   {
      // Slices only make sense in 3D
      Debug::StaticAssert< (DIMENSION == Dimensions::THREED) >();

      // Check size
      assert(k < this->mField.size());

      return this->mField.at(k);
   }

   template <typename TData, Dimensions::Type DIMENSION> inline typename SlicedLayout<TData, DIMENSION>::StorageType& SlicedLayout<TData, DIMENSION>::rSlice(const int k)
   {
      // Slices only make sense in 3D
      Debug::StaticAssert< (DIMENSION == Dimensions::THREED) >();

      // Check size
      assert(k < this->mField.size());

      return this->mField.at(k);
   }

   template <typename TData, Dimensions::Type DIMENSION> inline typename SlicedLayout<TData, DIMENSION>::StorageType::ConstColXpr SlicedLayout<TData, DIMENSION>::profile(const int j, const int k) const
   {
      // Profiles only make sense in 2D and 3D
      Debug::StaticAssert< (DIMENSION == Dimensions::TWOD || DIMENSION == Dimensions::THREED) >();

      // 3D version
      if(DIMENSION == Dimensions::THREED)
      {
         // Check size
         assert(k < this->mField.size());
         assert(j < this->mField.at(k).cols());

         return this->mField.at(k).col(j);

      // 2D version
      } else if(DIMENSION == Dimensions::TWOD)
      {
         // Check size
         assert(k == 0);
         assert(j < this->mField.size());

         return this->mField.at(j).col(0);
      }
   }

   template <typename TData, Dimensions::Type DIMENSION> inline typename SlicedLayout<TData, DIMENSION>::StorageType::ColXpr SlicedLayout<TData, DIMENSION>::rProfile(const int j, const int k)
   {
      // Profiles only make sense in 2D and 3D
      Debug::StaticAssert< (DIMENSION == Dimensions::TWOD || DIMENSION == Dimensions::THREED) >();

      // 3D version
      if(DIMENSION == Dimensions::THREED)
      {
         // Check size
         assert(k < this->mField.size());
         assert(j < this->mField.at(k).cols());

         return this->mField.at(k).col(j);

      // 2D version
      } else if(DIMENSION == Dimensions::TWOD)
      {
         // Check size
         assert(k == 0);
         assert(j < this->mField.size());

         return this->mField.at(j).col(0);
      }
   }

   template <typename TData, Dimensions::Type DIMENSION> inline const TData& SlicedLayout<TData, DIMENSION>::point(const int i, const int j, const int k) const
   {
      // 3D version
      if(DIMENSION == Dimensions::THREED)
      {
         // Check size
         assert(k < this->mField.size());
         assert(j < this->mField.at(k).cols());
         assert(i < this->mField.at(k).rows());

         return this->mField.at(k)(i,j);

      // 2D version
      } else if(DIMENSION == Dimensions::TWOD)
      {
         // Check size
         assert(k == 0);
         assert(j < this->mField.size());
         assert(i < this->mField.at(j).rows());

         return this->mField.at(j)(i,0);

      // 1D version
      } else if(DIMENSION == Dimensions::ONED)
      {
         // Check size
         assert(k == 0);
         assert(j == 0);
         assert(i < this->mField.at(0).rows());

         return this->mField.at(0)(i,0);
      }
   }

   template <typename TData, Dimensions::Type DIMENSION> inline TData& SlicedLayout<TData, DIMENSION>::rPoint(const int i, const int j, const int k)
   {
      // 3D version
      if(DIMENSION == Dimensions::THREED)
      {
         // Check size
         assert(k < this->mField.size());
         assert(j < this->mField.at(k).cols());
         assert(i < this->mField.at(k).rows());

         return this->mField.at(k)(i,j);

      // 2D version
      } else if(DIMENSION == Dimensions::TWOD)
      {
         // Check size
         assert(k == 0);
         assert(j < this->mField.size());
         assert(i < this->mField.at(j).rows());

         return this->mField.at(j)(i,0);

      // 1D version
      } else if(DIMENSION == Dimensions::ONED)
      {
         // Check size
         assert(k == 0);
         assert(j == 0);
         assert(i < this->mField.at(0).rows());

         return this->mField.at(0)(i,0);
      }
   }

   template <typename TData, Dimensions::Type DIMENSION> SlicedLayout<TData, DIMENSION>::SlicedLayout()
   {
   }

   template <typename TData, Dimensions::Type DIMENSION> void SlicedLayout<TData, DIMENSION>::initStorage(const ArrayI& dim1D, const ArrayI& dim2D, const int dim3D)
   {
      // Check for correct sizes of input data
      if(dim1D.size() != std::max(1, (dim3D > 0 ? dim3D : dim2D(0))) || dim2D.size() != std::max(1,dim3D))
      {
         throw Exception("SlicedLayout::initStorage", "Given resolution input data does not have the right size!");
      }

      // Check that initialisation is possible
      if((dim1D.array() < 0).any() || (dim2D.array() < 0).any() || dim3D < 0 || (dim2D.sum() == 0 && dim3D > 0))
      {
         throw Exception("SlicedLayout::initStorage", "Resolution information is not compatible with the sliced layout!");
      } else
      {
         // Check if it is 3D
         if(dim3D > 0)
         {
            // Initialise the storage
            this->mField.reserve(dim3D);
            for(int k = 0; k < dim3D; ++k)
            {
               this->mField.push_back(StorageType(dim1D(k), dim2D(k)));
            }

         // Check if it is 2D
         } else if(dim3D == 0 && dim2D.size() == 1 && dim2D(0) > 0)
         {
            // Initialise the storage for each set
            this->mField.reserve(dim2D(0));
            for(int j = 0; j < dim2D(0); ++j)
            {
               this->mField.push_back(StorageType(dim1D(j), 1));
            }

         // Check if it is 1D
         } else if(dim3D == 0 && dim2D.size() == 1 && dim2D(0) == 0 && dim1D.size() == 1 && dim1D(0) > 0)
         {
            // Initialise the storage Map
            this->mField.reserve(1);
            this->mField.push_back(StorageType(dim1D(0), 1));

         // Setup is not compatible
         } else
         {
            throw Exception("SlicedLayout::initStorage", "Could not determine the dimension of the scalar!");
         }
      }

      // Make sure everything is initialised to zero
      this->initialiseZeros();
   }

   template <typename TData, Dimensions::Type DIMENSION> void SlicedLayout<TData, DIMENSION>::initialiseZeros()
   {
      typename std::vector<typename SlicedLayout<TData, DIMENSION>::StorageType>::iterator it;
      for(it = this->mField.begin(); it !=this->mField.end(); it++)
      {
         it->setConstant(0.0);
      }
   }

   template <typename TData, Dimensions::Type DIMENSION> void SlicedLayout<TData, DIMENSION>::rescale(const MHDFloat scale)
   {
      typename std::vector<typename SlicedLayout<TData, DIMENSION>::StorageType>::iterator it;
      for(it = this->mField.begin(); it !=this->mField.end(); it++)
      {
         (*it) *= scale;
      }
   }

#ifdef GEOMHDISCC_STORAGEPROFILE
   template <typename TData, Dimensions::Type DIMENSION> MHDFloat SlicedLayout<TData, DIMENSION>::requiredStorage() const
   {
      MHDFloat mem = 0.0;

      typename std::vector<typename SlicedLayout<TData, DIMENSION>::StorageType>::const_iterator it;
      for(it = this->mField.begin(); it !=this->mField.end(); it++)
      {
         mem += it->size();
      }

      return static_cast<MHDFloat>(Debug::MemorySize<TData>::BYTES)*mem;
   }
#endif // GEOMHDISCC_STORAGEPROFILE

}
}

#endif // SLICEDLAYOUT_HPP
