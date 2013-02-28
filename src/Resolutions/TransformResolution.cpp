/** \file TransformResolution.cpp
 *  \brief Source of the resolution object for a transform
 */

// System includes
//
#include <assert.h>

// External includes
//

// Class include
//
#include "Resolutions/TransformResolution.hpp"

// Project includes
//

namespace GeoMHDiSCC {

   TransformResolution::TransformResolution(const std::vector<ArrayI>& fwd, const std::vector<ArrayI>& bwd, const std::vector<ArrayI>& idx2D, const ArrayI& idx3D)
      : mFwd(fwd), mBwd(bwd), mIdx2D(idx2D), mIdx3D(idx3D)
   {
   }

   TransformResolution::~TransformResolution()
   {
   }

   template <> int TransformResolution::dim<Dimensions::Data::DATF1D>(const int k) const
   {
      // Check for correct sizes
      assert(k >= 0);
      assert(static_cast<size_t>(k) < this->mFwd.size());

      return this->mFwd.at(k).rows();
   }

   template <> int TransformResolution::dim<Dimensions::Data::DATF1D>() const
   {
      return this->dim<Dimensions::Data::DATF1D>(0);
   }

   template <> int TransformResolution::dim<Dimensions::Data::DATB1D>(const int k) const
   {
      // Check for correct sizes
      assert(k >= 0);
      assert(static_cast<size_t>(k) < this->mBwd.size());

      return this->mBwd.at(k).rows();
   }

   template <> int TransformResolution::dim<Dimensions::Data::DATB1D>() const
   {
      return this->dim<Dimensions::Data::DATB1D>(0);
   }

   template <> int TransformResolution::dim<Dimensions::Data::DAT2D>(const int k) const
   {
      // Check for correct sizes
      assert(k >= 0);
      assert(static_cast<size_t>(k) < this->mIdx2D.size());

      return this->mIdx2D.at(k).size();
   }

   template <> int TransformResolution::dim<Dimensions::Data::DAT2D>() const
   {
      return this->dim<Dimensions::Data::DAT2D>(0);
   }

   template <> int TransformResolution::dim<Dimensions::Data::DAT3D>() const
   {
      return this->mIdx3D.size();
   }

   template <> int TransformResolution::idx<Dimensions::Data::DATF1D>(const int i, const int k) const
   {
      // Check for correct sizes
      assert(k >= 0);
      assert(static_cast<size_t>(k) < this->mFwd.size());
      assert(i < this->mFwd.at(k).size());

      return this->mFwd.at(k)(i);
   }

   template <> int TransformResolution::idx<Dimensions::Data::DATF1D>(const int i) const
   {
      return this->idx<Dimensions::Data::DATF1D>(i, 0);
   }

   template <> int TransformResolution::idx<Dimensions::Data::DATB1D>(const int i, const int k) const
   {
      // Check for correct sizes
      assert(k >= 0);
      assert(static_cast<size_t>(k) < this->mBwd.size());
      assert(i < this->mBwd.at(k).size());

      return this->mBwd.at(k)(i);
   }

   template <> int TransformResolution::idx<Dimensions::Data::DATB1D>(const int i) const
   {
      return this->idx<Dimensions::Data::DATB1D>(i, 0);
   }

   template <> int TransformResolution::idx<Dimensions::Data::DAT2D>(const int i, const int k) const
   {
      // Check for correct sizes
      assert(k >= 0);
      assert(static_cast<size_t>(k) < this->mIdx2D.size());
      assert(i < this->mIdx2D.at(k).size());

      return this->mIdx2D.at(k)(i);
   }

   template <> int TransformResolution::idx<Dimensions::Data::DAT2D>(const int i) const
   {
      return this->idx<Dimensions::Data::DAT2D>(i, 0);
   }

   template <> int TransformResolution::idx<Dimensions::Data::DAT3D>(const int i) const
   {
      // Check for correct sizes
      assert(i < this->mIdx3D.size());

      return this->mIdx3D(i);
   }

   Datatypes::SharedScalarFieldSetupType TransformResolution::spFwdSetup() const
   {
      // Iterator for the vector based storages
      std::vector<ArrayI>::const_iterator   vIt;
      
      // Get forward dimensions
      SharedArrayI   spDim1D(new ArrayI(this->dim<Dimensions::Data::DAT3D>()));
      int i = 0;
      for(vIt = this->mFwd.begin(); vIt != this->mFwd.end(); ++vIt, ++i)
      {
         (*spDim1D)(i) = vIt->size();
      }

      // Get 2D dimensions
      SharedArrayI   spDim2D(new ArrayI(this->dim<Dimensions::Data::DAT3D>()));
      i = 0;
      for(vIt = this->mIdx2D.begin(); vIt != this->mIdx2D.end(); ++vIt, ++i)
      {
         (*spDim2D)(i) = vIt->size();
      }

      return Datatypes::SharedScalarFieldSetupType(new Datatypes::ScalarFieldSetupType(spDim1D, spDim2D, this->dim<Dimensions::Data::DAT3D>()));
   }

   Datatypes::SharedScalarFieldSetupType TransformResolution::spBwdSetup() const
   {
      // Iterator for the vector based storages
      std::vector<ArrayI>::const_iterator   vIt;
      
      // Get backward dimensions
      SharedArrayI   spDim1D(new ArrayI(this->dim<Dimensions::Data::DAT3D>()));
      int i = 0;
      for(vIt = this->mBwd.begin(); vIt != this->mBwd.end(); ++vIt, ++i)
      {
         (*spDim1D)(i) = vIt->size();
      }

      // Get 2D dimensions
      SharedArrayI   spDim2D(new ArrayI(this->dim<Dimensions::Data::DAT3D>()));
      i = 0;
      for(vIt = this->mIdx2D.begin(); vIt != this->mIdx2D.end(); ++vIt, ++i)
      {
         (*spDim2D)(i) = vIt->size();
      }

      return Datatypes::SharedScalarFieldSetupType(new Datatypes::ScalarFieldSetupType(spDim1D, spDim2D, this->dim<Dimensions::Data::DAT3D>()));
   }

}
