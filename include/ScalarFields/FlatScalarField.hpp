/** \file FlatScalarField.hpp
 *  \brief Base for a  scalar field with flat data layout
 *
 *  \mhdBug Needs test
 */

#ifndef FLATSCALARFIELD_HPP
#define FLATSCALARFIELD_HPP

// Debug includes
//
#include "StorageProfiler/StorageProfilerMacro.h"
#include "StaticAsserts/StaticAssert.hpp"
#include "Exceptions/Exception.hpp"

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
    * @brief Base for a  scalar field with flat data layout
    */
   template <typename TData, Dimensions::Type DIMENSION> class FlatScalarField
   {
      public:
         /// Dimension of the scalar field
         static const Dimensions::Type FieldDimension = DIMENSION;

         /// Typedef for the coefficient type
         typedef ScalarFieldSetup<DIMENSION> SetupType;

         /// Typedef for the coefficient type
         typedef SharedPtrMacro<ScalarFieldSetup<DIMENSION> > SharedSetupType;

         /// Typedef for the coefficient type
         typedef TData PointType;

         /// Typedef for the storage type
         typedef typename Eigen::Matrix<PointType, Eigen::Dynamic, Eigen::Dynamic> StorageType;

         /// Typedef for the profiles of the storage type
         typedef typename StorageType::ColXpr  ProfileType;

         /// Typedef for the slices of the storage type
         typedef typename Eigen::Block<StorageType>  SliceType;

         /**
          * @brief Constructor
          */
         explicit FlatScalarField(SharedPtrMacro<ScalarFieldSetup<DIMENSION> > spSetup);

         /**
          * @brief Destructor
          */
         ~FlatScalarField();

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
          * @brief Add to 2D slice of the field
          *
          * @param sl   Slice values
          * @param k    Index of the slice
          */
         template <typename Derived> void addSlice(const Eigen::MatrixBase<Derived>& sl, const int k);

         /**
          * @brief Substract from 2D slice of the field
          *
          * @param sl   Slice values
          * @param k    Index of the slice
          */
         template <typename Derived> void subSlice(const Eigen::MatrixBase<Derived>& sl, const int k);

         /**
          * @brief Get internal storage field data
          */
         const StorageType& data() const;

         /**
          * @brief Set internal storage field data
          */
         template <typename Derived> void setData(const Eigen::MatrixBase<Derived>& field);

         /**
          * @brief Add to internal storage field data
          */
         template <typename Derived> void addData(const Eigen::MatrixBase<Derived>& field);

         /**
          * @brief Substract from internal storage field data
          */
         template <typename Derived> void subData(const Eigen::MatrixBase<Derived>& field);

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

         /**
          * @brief Set internal storage field data
          *
          * \warning This routine should only be used in exceptional cases. Use setData, addData, subData when you can!
          */
         StorageType& rData();

         /**
          * @brief Set internal storage point data
          *
          * \warning This routine should only be used in exceptional cases. Use setPoint!
          */
         PointType& rPoint(const int i, const int j = 0, const int k = 0);

      protected:

      private:
         /**
          * @brief Setup object for the scalar field
          */
         SharedPtrMacro<ScalarFieldSetup<DIMENSION> > mspSetup;

         /**
          * @brief Field values in shared flat storage
          */
         SharedPtrMacro<StorageType>  mspField;
   };

   template <typename TData, Dimensions::Type DIMENSION> inline typename FlatScalarField<TData, DIMENSION>::PointType FlatScalarField<TData, DIMENSION>::point(const int i, const int j, const int k) const
   {
      return (*this->mspField)(i,this->mspSetup->colIdx(j,k));
   }

   template <typename TData, Dimensions::Type DIMENSION> inline typename FlatScalarField<TData, DIMENSION>::ProfileType FlatScalarField<TData, DIMENSION>::profile(const int j, const int k) const
   {
      // Profiles only make sense in 2D and 3D
      Debug::StaticAssert< (DIMENSION == Dimensions::TWOD) || (DIMENSION == Dimensions::THREED) >();

      return this->mspField->col(this->mspSetup->colIdx(j,k));
   }

   template <typename TData, Dimensions::Type DIMENSION> inline typename FlatScalarField<TData, DIMENSION>::SliceType FlatScalarField<TData, DIMENSION>::slice(const int k) const
   {
      // Slices only make sense in 3D
      Debug::StaticAssert< (DIMENSION == Dimensions::THREED) >();

      return this->mspField->block(0,this->mspSetup->blockIdx(k), this->mspSetup->blockRows(k), this->mspSetup->blockCols(k));
   }

   template <typename TData, Dimensions::Type DIMENSION> inline const typename FlatScalarField<TData, DIMENSION>::StorageType& FlatScalarField<TData, DIMENSION>::data() const
   {
      return *this->mspField;
   }

   template <typename TData, Dimensions::Type DIMENSION> inline typename FlatScalarField<TData, DIMENSION>::StorageType& FlatScalarField<TData, DIMENSION>::rData()
   {
      return *this->mspField;
   }

   template <typename TData, Dimensions::Type DIMENSION> void FlatScalarField<TData, DIMENSION>::setPoint(const FlatScalarField<TData, DIMENSION>::PointType pt, const int i, const int j, const int k)
   {
      (*this->mspField)(i,this->mspSetup->colIdx(j,k)) = pt;
   }

   template <typename TData, Dimensions::Type DIMENSION> inline typename FlatScalarField<TData, DIMENSION>::PointType& FlatScalarField<TData, DIMENSION>::rPoint(const int i, const int j, const int k)
   {
      return (*this->mspField)(i,this->mspSetup->colIdx(j,k));
   }

   template <typename TData, Dimensions::Type DIMENSION> template<typename Derived> void FlatScalarField<TData, DIMENSION>::setProfile(const Eigen::MatrixBase<Derived>& pf, const int j, const int k)
   {
      // Profiles only make sense in 2D and 3D
      Debug::StaticAssert< (DIMENSION == Dimensions::TWOD) || (DIMENSION == Dimensions::THREED) >();

      this->mspField->col(this->mspSetup->colIdx(j,k)) = pf;
   }

   template <typename TData, Dimensions::Type DIMENSION> template<typename Derived> void FlatScalarField<TData, DIMENSION>::setSlice(const Eigen::MatrixBase<Derived>& sl, const int k)
   {
      // Slices only make sense in 3D
      Debug::StaticAssert< (DIMENSION == Dimensions::THREED) >();

      this->mspField->block(0, this->mspSetup->blockIdx(k), this->mspSetup->blockRows(k), this->mspSetup->blockCols(k)) = sl;
   }

   template <typename TData, Dimensions::Type DIMENSION> template<typename Derived> void FlatScalarField<TData, DIMENSION>::addSlice(const Eigen::MatrixBase<Derived>& sl, const int k)
   {
      // Slices only make sense in 3D
      Debug::StaticAssert< (DIMENSION == Dimensions::THREED) >();

      this->mspField->block(0, this->mspSetup->blockIdx(k), this->mspSetup->blockRows(k), this->mspSetup->blockCols(k)) += sl;
   }

   template <typename TData, Dimensions::Type DIMENSION> template<typename Derived> void FlatScalarField<TData, DIMENSION>::subSlice(const Eigen::MatrixBase<Derived>& sl, const int k)
   {
      // Slices only make sense in 3D
      Debug::StaticAssert< (DIMENSION == Dimensions::THREED) >();

      this->mspField->block(0, this->mspSetup->blockIdx(k), this->mspSetup->blockRows(k), this->mspSetup->blockCols(k)) -= sl;
   }

   template <typename TData, Dimensions::Type DIMENSION> template<typename Derived> void FlatScalarField<TData, DIMENSION>::setData(const Eigen::MatrixBase<Derived>& field)
   {
      this->mspField->block(0, 0, this->mspSetup->dataRows(), this->mspSetup->dataCols()) = field;
   }

   template <typename TData, Dimensions::Type DIMENSION> template<typename Derived> void FlatScalarField<TData, DIMENSION>::addData(const Eigen::MatrixBase<Derived>& field)
   {
      this->mspField->block(0, 0, this->mspSetup->dataRows(), this->mspSetup->dataCols()) += field;
   }

   template <typename TData, Dimensions::Type DIMENSION> template<typename Derived> void FlatScalarField<TData, DIMENSION>::subData(const Eigen::MatrixBase<Derived>& field)
   {
      this->mspField->block(0, 0, this->mspSetup->dataRows(), this->mspSetup->dataCols()) -= field;
   }

   template <typename TData, Dimensions::Type DIMENSION> FlatScalarField<TData, DIMENSION>::FlatScalarField(SharedPtrMacro<ScalarFieldSetup<DIMENSION> > spSetup)
      : mspSetup(spSetup), mspField(new StorageType(spSetup->dataRows(), spSetup->dataCols()))
   {
   }

   template <typename TData, Dimensions::Type DIMENSION> FlatScalarField<TData, DIMENSION>::~FlatScalarField()
   {
   }

   template <typename TData, Dimensions::Type DIMENSION> void FlatScalarField<TData, DIMENSION>::setZeros()
   {
      this->mspField->setConstant(0.0);
   }

   template <typename TData, Dimensions::Type DIMENSION> void FlatScalarField<TData, DIMENSION>::rescale(const MHDFloat scale)
   {
      *this->mspField *= scale;
   }

#ifdef GEOMHDISCC_STORAGEPROFILE
   template <typename TData, Dimensions::Type DIMENSION> MHDFloat FlatScalarField<TData, DIMENSION>::requiredStorage() const
   {
      MHDFloat mem = this->mspField->size();

      return static_cast<MHDFloat>(Debug::MemorySize<TData>::BYTES)*mem;
   }
#endif // GEOMHDISCC_STORAGEPROFILE

}
}

#endif // FLATSCALARFIELD_HPP
