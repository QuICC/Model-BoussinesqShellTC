/** \file FlatLayout.hpp
 *  \brief Base for a  scalar field with flat data layout
 */

#ifndef FLATLAYOUT_HPP
#define FLATLAYOUT_HPP

// Configuration includes
//
#include "StorageProfiler/StorageProfilerMacro.h"

// System includes
//
#include <assert.h>
#include <vector>
#include <algorithm>

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Exceptions/Exception.hpp"
#include "Enums/Dimensions.hpp"
#include "StaticAssert/StaticAssert.hpp"

namespace GeoMHDiSCC {

namespace Datatypes {

   /**
    * @brief Base for a  scalar field with flat data layout
    */
   template <typename TData, Dimensions::Type DIMENSION> class FlatLayout
   {
      public:
         /// Typedef for the coefficient type
         typedef TData CoefficientType;

         /// Typedef for the storage type
         typedef Eigen::Matrix<TData, Eigen::Dynamic, Eigen::Dynamic> StorageType;

         /// Typedef for the Map to the storage type
         typedef Eigen::Map<StorageType>  MapStorageType;

         /**
          * @brief Constructor
          */
         FlatLayout();

         /**
          * @brief Copy constructor is required because of the vector of Map objects
          *
          * @param oldLayout Layout to copy over
          */
         FlatLayout(const FlatLayout& oldLayout);

         /**
          * @brief Destructor
          */
         virtual ~FlatLayout();

         /**
          * @brief Get sliced field data
          */
         const std::vector<MapStorageType>& sliced() const;

         /**
          * @brief Set sliced field data
          */
         std::vector<MapStorageType>& rSliced();

         /**
          * @brief Get internal storage field data
          */
         const StorageType& data() const;

         /**
          * @brief Set internal storage field data
          */
         StorageType& rData();

         /**
          * @brief Get full field data as a single column
          */
         const MapStorageType& linear() const;

         /**
          * @brief Set full field data as a single column
          */
         MapStorageType& rLinear();

         /**
          * @brief Get a 2D slice of the field
          *
          * @param k Index of the slice
          */
         const MapStorageType& slice(const int k) const;

         /**
          * @brief Set a 2D slice of the field
          *
          * @param k Index of the slice
          */
         MapStorageType& rSlice(const int k);

         /**
          * @brief Get a 1D profiled of the field
          *
          * @param j Index of the profile
          * @param k Index of the slice
          */
         typename MapStorageType::ConstColXpr profile(const int j, const int k = 0) const;

         /**
          * @brief Set a 1D profiled of the field
          *
          * @param j Index of the profile
          * @param k Index of the slice
          */
         typename MapStorageType::ColXpr rProfile(const int j, const int k = 0);

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

         /**
          * @brief Field values in flat storage
          */
         StorageType  mField;

         /**
          * @brief Linear map to the data storage
          */
         MapStorageType mLinearMap;

         /**
          * @brief 
         
         /**
          * @brief Create map of flat storage to sliced storage
          */
         std::vector<MapStorageType>  mFieldMap;

      protected:
         /**
          * @brief Initialise the storage
          * @param dim1D Size in the first dimension
          * @param dim2D Size in the second dimension
          * @param dim3D Size in the third dimension
          */
         void initStorage(const ArrayI& dim1D, const ArrayI& dim2D, const int dim3D);

      private:
   };

   template <typename TData, Dimensions::Type DIMENSION> inline const std::vector<typename FlatLayout<TData, DIMENSION>::MapStorageType>& FlatLayout<TData, DIMENSION>::sliced() const
   {
      return this->mFieldMap;
   }

   template <typename TData, Dimensions::Type DIMENSION> inline std::vector<typename FlatLayout<TData, DIMENSION>::MapStorageType>& FlatLayout<TData, DIMENSION>::rSliced()
   {
      return this->mFieldMap;
   }

   template <typename TData, Dimensions::Type DIMENSION> inline const typename FlatLayout<TData, DIMENSION>::StorageType& FlatLayout<TData, DIMENSION>::data() const
   {
      return this->mField;
   }

   template <typename TData, Dimensions::Type DIMENSION> inline typename FlatLayout<TData, DIMENSION>::StorageType& FlatLayout<TData, DIMENSION>::rData()
   {
      return this->mField;
   }

   template <typename TData, Dimensions::Type DIMENSION> inline const typename FlatLayout<TData, DIMENSION>::MapStorageType& FlatLayout<TData, DIMENSION>::linear() const
   {
      return this->mLinearMap;
   }

   template <typename TData, Dimensions::Type DIMENSION> inline typename FlatLayout<TData, DIMENSION>::MapStorageType& FlatLayout<TData, DIMENSION>::rLinear()
   {
      return this->mLinearMap;
   }

   template <typename TData, Dimensions::Type DIMENSION> inline const typename FlatLayout<TData, DIMENSION>::MapStorageType& FlatLayout<TData, DIMENSION>::slice(const int k) const
   {
      // Slices only make sense in 3D
      StaticAssert< (DIMENSION == Dimensions::THREED) >();

      // Check size
      assert(k < this->mFieldMap.size());

      return this->mFieldMap.at(k);
   }

   template <typename TData, Dimensions::Type DIMENSION> inline typename FlatLayout<TData, DIMENSION>::MapStorageType& FlatLayout<TData, DIMENSION>::rSlice(const int k)
   {
      // Slices only make sense in 3D
      StaticAssert< (DIMENSION == Dimensions::THREED) >();

      // Check size
      assert(k < this->mFieldMap.size());

      return this->mFieldMap.at(k);
   }

   template <typename TData, Dimensions::Type DIMENSION> inline typename FlatLayout<TData, DIMENSION>::MapStorageType::ConstColXpr FlatLayout<TData, DIMENSION>::profile(const int j, const int k) const
   {
      // Profiles only make sense in 2D and 3D
      StaticAssert< (DIMENSION == Dimensions::TWOD || DIMENSION == Dimensions::THREED) >();

      // 3D version
      if(DIMENSION == Dimensions::THREED)
      {
         // Check size
         assert(k < this->mFieldMap.size());
         assert(j < this->mFieldMap.at(k).cols());

         return this->mFieldMap.at(k).col(j);

      // 2D version
      } else if(DIMENSION == Dimensions::TWOD)
      {
         // Check size
         assert(k == 0);
         assert(j < this->mFieldMap.size());

         return this->mFieldMap.at(j).col(0);
      }
   }

   template <typename TData, Dimensions::Type DIMENSION> inline typename FlatLayout<TData, DIMENSION>::MapStorageType::ColXpr FlatLayout<TData, DIMENSION>::rProfile(const int j, const int k)
   {
      // Profiles only make sense in 2D and 3D
      StaticAssert< (DIMENSION == Dimensions::TWOD || DIMENSION == Dimensions::THREED) >();

      // 3D version
      if(DIMENSION == Dimensions::THREED)
      {
         // Check size
         assert(k < this->mFieldMap.size());
         assert(j < this->mFieldMap.at(k).cols());

         return this->mFieldMap.at(k).col(j);

      // 2D version
      } else if(DIMENSION == Dimensions::TWOD)
      {
         // Check size
         assert(k == 0);
         assert(j < this->mFieldMap.size());

         return this->mFieldMap.at(j).col(0);
      }
   }

   template <typename TData, Dimensions::Type DIMENSION> inline const TData& FlatLayout<TData, DIMENSION>::point(const int i, const int j, const int k) const
   {
      // 3D version
      if(DIMENSION == Dimensions::THREED)
      {
         // Check size
         assert(k < this->mFieldMap.size());
         assert(j < this->mFieldMap.at(k).cols());
         assert(i < this->mFieldMap.at(k).rows());

         return this->mFieldMap.at(k)(i,j);

      // 2D version
      } else if(DIMENSION == Dimensions::TWOD)
      {
         // Check size
         assert(k == 0);
         assert(j < this->mFieldMap.size());
         assert(i < this->mFieldMap.at(j).rows());

         return this->mFieldMap.at(j)(i,0);

      // 1D version
      } else if(DIMENSION == Dimensions::ONED)
      {
         // Check size
         assert(k == 0);
         assert(j == 0);
         assert(i < this->mFieldMap.at(0).rows());

         return this->mFieldMap.at(0)(i,0);
      }

   }

   template <typename TData, Dimensions::Type DIMENSION> inline TData& FlatLayout<TData, DIMENSION>::rPoint(const int i, const int j, const int k)
   {
      // 3D version
      if(DIMENSION == Dimensions::THREED)
      {
         // Check size
         assert(k < this->mFieldMap.size());
         assert(j < this->mFieldMap.at(k).cols());
         assert(i < this->mFieldMap.at(k).rows());

         return this->mFieldMap.at(k)(i,j);

      // 2D version
      } else if(DIMENSION == Dimensions::TWOD)
      {
         // Check size
         assert(k == 0);
         assert(j < this->mFieldMap.size());
         assert(i < this->mFieldMap.at(j).rows());

         return this->mFieldMap.at(j)(i,0);

      // 1D version
      } else if(DIMENSION == Dimensions::ONED)
      {
         // Check size
         assert(k == 0);
         assert(j == 0);
         assert(i < this->mFieldMap.at(0).rows());

         return this->mFieldMap.at(0)(i,0);
      }
   }

   template <typename TData, Dimensions::Type DIMENSION> FlatLayout<TData, DIMENSION>::FlatLayout()
      : mField(1,1), mLinearMap(0,0,0)
   {
   }

   template <typename TData, Dimensions::Type DIMENSION> FlatLayout<TData, DIMENSION>::FlatLayout(const FlatLayout<TData, DIMENSION>& oldLayout)
      : mField(oldLayout.mField), mLinearMap(0,0,0)
   {
      // Initialise the map vector
      this->mFieldMap.resize(oldLayout.mFieldMap.size(), MapStorageType(0,0,0));

      // Put in the right data pointer 
      int shift = 0;
      for(int k = 0; k < this->mFieldMap.size(); ++k)
      {
         new (&this->mFieldMap[k]) MapStorageType(this->mField.data() + shift, oldLayout.mFieldMap.at(k).rows(), oldLayout.mFieldMap.at(k).cols());
         shift += oldLayout.mFieldMap.at(k).size();
      }

      // Recreate the linear map
      new (&this->mLinearMap) MapStorageType(this->mField.data(), oldLayout.mField.size(), 1);
   }

   template <typename TData, Dimensions::Type DIMENSION> FlatLayout<TData, DIMENSION>::~FlatLayout()
   {
   }

   template <typename TData, Dimensions::Type DIMENSION> void FlatLayout<TData, DIMENSION>::initStorage(const ArrayI& dim1D, const ArrayI& dim2D, const int dim3D)
   {
      // Check for correct sizes of input data
      if(dim1D.size() != std::max(1,(dim3D > 0 ? dim3D : dim2D(0))) || dim2D.size() != std::max(1,dim3D))
      {
         throw Exception("FlatLayout::initStorage", "Given resolution input data does not have the right size!");
      }

      // Check that initialisation is possible
      if((dim1D.array() < 0).any() || (dim2D.array() < 0).any() || dim3D < 0 || (dim2D.sum() == 0 && dim3D > 0) || !(dim1D.array() == dim1D(0)).all())
      {
         throw Exception("FlatLayout::initStorage", "Resolution information is not compatible with the flat layout!");
      } else
      {
         // Check if it is 3D
         if(dim3D > 0)
         {
            // Initialise the flat storage matrix
            this->mField.resize(dim1D(0), dim2D.sum());
            this->mField.setConstant(0.0);

            // Create the linear mapping
            new (&this->mLinearMap) MapStorageType(this->mField.data(), this->mField.size(), 1);

            // Initialise the map vector
            this->mFieldMap.resize(dim3D, MapStorageType(0,0,0));

            // Put in the right data pointer 
            int shift = 0;
            for(int k = 0; k < dim3D; ++k)
            {
               new (&this->mFieldMap[k]) MapStorageType(this->mField.data() + shift, dim1D(0), dim2D(k));
               shift += dim1D(0)*dim2D(k);
            }

         // Check if it is 2D
         } else if(dim3D == 0 && dim2D.size() == 1 && dim2D(0) > 0)
         {
            // Initialise the flat storage matrix
            this->mField.resize(dim1D(0), dim2D(0));
            this->mField.setConstant(0.0);

            // Create the linear mapping
            new (&this->mLinearMap) MapStorageType(this->mField.data(), this->mField.size(), 1);

            // Initialise the map vector
            this->mFieldMap.resize(dim2D(0), MapStorageType(0,0,0));

            // Put in the right data pointer
            int shift = 0;
            for(int j = 0; j < dim2D(0); ++j)
            {
               new (&this->mFieldMap[j]) MapStorageType(this->mField.data() + shift, dim1D(j), 1);
               shift += dim1D(j);
            }

         // Check if it is 1D
         } else if(dim3D == 0 && dim2D.size() == 1 && dim2D(0) == 0 && dim1D.size() == 1 && dim1D(0) > 0)
         {
            // Initialise the flat storage matrix
            this->mField.resize(dim1D(0), 1);
            this->mField.setConstant(0.0);

            // Create the linear mapping
            new (&this->mLinearMap) MapStorageType(this->mField.data(), this->mField.size(), 1);

            // Initialise the map vector
            this->mFieldMap.resize(dim1D(0), MapStorageType(0,0,0));

            // Put in the right data pointer 
            new (&this->mFieldMap[0]) MapStorageType(this->mField.data(), dim1D(0), 1);

         // Setup is not compatible
         } else
         {
            throw Exception("FlatLayout::initStorage", "Could not determine the dimension of the scalar!");
         }
      }

      // Make sure everything is initialised to zero
      this->initialiseZeros();
   }

   template <typename TData, Dimensions::Type DIMENSION> void FlatLayout<TData, DIMENSION>::initialiseZeros()
   {
      this->mField.setConstant(0.0);
   }

   template <typename TData, Dimensions::Type DIMENSION> void FlatLayout<TData, DIMENSION>::rescale(const MHDFloat scale)
   {
      this->mField *= scale;
   }

#ifdef GEOMHDISCC_STORAGEPROFILE
   template <typename TData, Dimensions::Type DIMENSION> MHDFloat FlatLayout<TData, DIMENSION>::requiredStorage() const
   {
      MHDFloat mem = this->mField.size();

      return static_cast<MHDFloat>(MemorySize<TData>::BYTES)*mem;
   }
#endif // GEOMHDISCC_STORAGEPROFILE

}
}

#endif // FLATLAYOUT_HPP
