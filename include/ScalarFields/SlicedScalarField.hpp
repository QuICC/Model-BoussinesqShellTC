/** 
 * @file SlicedScalarField.hpp
 * @brief Base for a  scalar field with sliced data layout
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 * @version 0.9.0
 * @date 2013-09-11
 *
 *  \mhdBug Needs corrected implementation
 */

#ifndef SLICEDSCALARFIELD_HPP
#define SLICEDSCALARFIELD_HPP

// Debug includes
//
#include "StorageProfiler/StorageProfilerMacro.h"
#include "Exceptions/Exception.hpp"
#include "StaticAsserts/StaticAssert.hpp"

// Configuration includes
//
#include "SmartPointers/SharedPtrMacro.h"

// System includes
//
#include <cassert>
#include <vector>

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Enums/Dimensions.hpp"
#include "ScalarFields/ScalarFieldSetup.hpp"

namespace GeoMHDiSCC {

namespace Datatypes {

   /**
    * @brief Base for a  scalar field with sliced data layout
    */
   template <typename TData, Dimensions::Type DIMENSION> class SlicedScalarField
   {
      public:
         /// Typedef for the coefficient type
         typedef ScalarFieldSetup<DIMENSION> SetupType;

         /// Typedef for the coefficient type
         typedef TData PointType;

         /// Typedef for the storage type
         typedef typename std::vector<Eigen::Matrix<PointType, Eigen::Dynamic, Eigen::Dynamic> > StorageType;

         /// Typedef for the profiles of the storage type
         typedef typename Eigen::Matrix<PointType, Eigen::Dynamic, Eigen::Dynamic>::ColXpr  ProfileType;

         /// Typedef for the slices of the storage type
         typedef typename Eigen::Matrix<PointType, Eigen::Dynamic, Eigen::Dynamic>  SliceType;

         /**
          * @brief Constructor
          */
         explicit SlicedScalarField(SharedPtrMacro<ScalarFieldSetup<DIMENSION> > spSetup);

         /**
          * @brief Destructor
          */
         ~SlicedScalarField();

         /**
          * @brief Get a point value of the field
          *
          * @param i Index of the point
          * @param j Index of the profile
          * @param k Index of the slice
          */
         PointType point(const int i, const int j = 0, const int k = 0) const;

         /**
          * @brief Get a 1D profile of the field 
          *
          * A profile is defined as all value for fixed indexes in the second and third dimension.
          *
          * @param j Index of the profile
          * @param k Index of the slice
          */
         ProfileType profile(const int j, const int k = 0) const;

         /**
          * @brief Get a 2D slice of the field
          *
          * A slice is the matrix of values for a fixed index in the third dimension
          *
          * @param k Index of the slice
          */
         SliceType slice(const int k) const;

         /**
          * @brief Set a point value of the field
          *
          * @param i Index of the point
          * @param j Index of the profile
          * @param k Index of the slice
          */
         void setPoint(const PointType pt, const int i, const int j = 0, const int k = 0);

         /**
          * @brief Set a profile of the field
          *
          * @param pf   Profile values
          * @param j    Index of the profile
          * @param k    Index of the slice
          */
         template <typename Derived> void setProfile(const Eigen::MatrixBase<Derived>& pf, const int j, const int k = 0);

         /**
          * @brief Set a 2D slice of the field
          *
          * @param sl   Slice values
          * @param k    Index of the slice
          */
         template <typename Derived> void setSlice(const Eigen::MatrixBase<Derived>& sl, const int k);

         /**
          * @brief Get internal storage field data
          */
         const StorageType& data() const;

         /**
          * @brief Set internal storage field data
          */
         template <typename Derived> void setData(const Eigen::MatrixBase<Derived>& field);

         /**
          * @brief Set the complete field to zero
          */
         void setZeros();

         /**
          * @brief Rescale the complete field by a real coefficient
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

      private:
         /**
          * @brief Setup object for the scalar field
          */
         SharedPtrMacro<ScalarFieldSetup<DIMENSION> > mspSetup;

         /**
          * @brief Field values
          */
         SharedPtrMacro<StorageType>  mspField;
   };

   template <typename TData, Dimensions::Type DIMENSION> inline typename SlicedScalarField<TData, DIMENSION>::PointType SlicedScalarField<TData, DIMENSION>::point(const int i, const int j, const int k) const
   {
      return this->mspField->at(k)(i,j);
   }

   template <typename TData, Dimensions::Type DIMENSION> inline typename SlicedScalarField<TData, DIMENSION>::ProfileType SlicedScalarField<TData, DIMENSION>::profile(const int j, const int k) const
   {
      // Profiles only make sense in 2D and 3D
      Debug::StaticAssert< (DIMENSION == Dimensions::TWOD) || (DIMENSION == Dimensions::THREED) >();

      return this->mspField->at(k).col(j);
   }

   template <typename TData, Dimensions::Type DIMENSION> inline typename SlicedScalarField<TData, DIMENSION>::SliceType SlicedScalarField<TData, DIMENSION>::slice(const int k) const
   {
      // Slices only make sense in 3D
      Debug::StaticAssert< (DIMENSION == Dimensions::THREED) >();

      return this->mField->at(k);
   }

   template <typename TData, Dimensions::Type DIMENSION> inline const typename SlicedScalarField<TData, DIMENSION>::StorageType& SlicedScalarField<TData, DIMENSION>::data() const
   {
      return *this->mspField;
   }

   template <typename TData, Dimensions::Type DIMENSION> void SlicedScalarField<TData, DIMENSION>::setPoint(const SlicedScalarField<TData, DIMENSION>::PointType pt, const int i, const int j, const int k)
   {
      this->mspField->at(k)(i,j) = pt;
   }

   template <typename TData, Dimensions::Type DIMENSION> template<typename Derived> void SlicedScalarField<TData, DIMENSION>::setProfile(const Eigen::MatrixBase<Derived>& pf, const int j, const int k)
   {
      // Profiles only make sense in 2D and 3D
      Debug::StaticAssert< (DIMENSION == Dimensions::TWOD) || (DIMENSION == Dimensions::THREED) >();

      this->mspField->at(k).col(j) = pf;
   }

   template <typename TData, Dimensions::Type DIMENSION> template<typename Derived> void SlicedScalarField<TData, DIMENSION>::setSlice(const Eigen::MatrixBase<Derived>& sl, const int k)
   {
      // Slices only make sense in 3D
      Debug::StaticAssert< (DIMENSION == Dimensions::THREED) >();

      this->mspField->at(k) = sl;
   }

   template <typename TData, Dimensions::Type DIMENSION> template<typename Derived> void SlicedScalarField<TData, DIMENSION>::setData(const Eigen::MatrixBase<Derived>& field)
   {
      *this->mspField = field;
   }

   template <typename TData, Dimensions::Type DIMENSION> SlicedScalarField<TData, DIMENSION>::SlicedScalarField()
   {
   }

   template <typename TData, Dimensions::Type DIMENSION> SlicedScalarField<TData, DIMENSION>::~SlicedScalarField()
   {
   }

   template <typename TData, Dimensions::Type DIMENSION> void SlicedScalarField<TData, DIMENSION>::setZeros()
   {
      // Get iterator
      StorageType::iterator   it;

      // Iterate over all slices
      for(it = this->mspField->begin(); it != this->mspField->end(); ++it)
      {
         it->setConstant(0.0);
      }
   }

   template <typename TData, Dimensions::Type DIMENSION> void SlicedScalarField<TData, DIMENSION>::rescale(const MHDFloat scale)
   {
      // Get iterator
      StorageType::iterator   it;

      // Iterate over all slices
      for(it = this->mspField->begin(); it != this->mspField->end(); ++it)
      {
         (*it) *= scale;
      }
   }

#ifdef GEOMHDISCC_STORAGEPROFILE
   template <typename TData, Dimensions::Type DIMENSION> MHDFloat SlicedScalarField<TData, DIMENSION>::requiredStorage() const
   {
      MHDFloat mem = 0.0;

      StorageType::const_iterator it;
      for(it = this->mspField->begin(); it != this->mspField->end(); ++it)
      {
         mem += it->size();
      }

      return static_cast<MHDFloat>(Debug::MemorySize<TData>::BYTES)*mem;
   }
#endif // GEOMHDISCC_STORAGEPROFILE

}
}

#endif // SLICEDSCALARFIELD_HPP
