/** 
 * @file TransformResolution.cpp
 * @brief Source of the resolution object for a transform
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// System includes
//
#include <cassert>

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
      assert(this->mFwd.size() > 0);
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
      assert(this->mBwd.size() > 0);
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
      assert(this->mIdx2D.size() > 0);
      assert(static_cast<size_t>(k) < this->mIdx2D.size());

      return this->mIdx2D.at(k).size();
   }

   template <> int TransformResolution::dim<Dimensions::Data::DAT2D>() const
   {
      return this->dim<Dimensions::Data::DAT2D>(0);
   }

   template <> int TransformResolution::dim<Dimensions::Data::DAT3D>() const
   {
      // Check for correct size
      assert(this->mIdx3D.size() > 0);

      return this->mIdx3D.size();
   }

   template <> int TransformResolution::idx<Dimensions::Data::DATF1D>(const int i, const int k) const
   {
      // Check for correct sizes
      assert(k >= 0);
      assert(static_cast<size_t>(k) < this->mFwd.size());
      assert(this->mFwd.at(k).size() > 0);
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
      assert(this->mBwd.at(k).size() > 0);
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
      assert(this->mIdx2D.at(k).size() > 0);
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
      assert(this->mIdx3D.size() > 0);
      assert(i < this->mIdx3D.size());

      return this->mIdx3D(i);
   }

   ArrayI TransformResolution::mode(const int i) const
   {
      ArrayI mode(4); // 0 -> index 3D, 1 -> index 2D, 2 -> mode 3D, 3 -> mode 2D

      int current = 0;
      for(mode(0) = 0; mode(0) < this->mIdx3D.size(); ++mode(0))
      {
         if(i - current < this->mIdx2D.at(mode(0)).size())
         {
            mode(1) = i-current;
            mode(2) = this->mIdx3D(mode(0));
            mode(3) = this->mIdx2D.at(mode(0))(mode(1));
            current = i;
            break;
         } else
         {
            current += this->mIdx2D.at(mode(0)).size();
         }
      }

      // Safety assert
      assert(current == i);

      return mode;
   }

}
