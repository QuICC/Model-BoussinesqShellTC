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

   int TransformResolution::dimFwd(const int j, const int k) const
   {
      if(this->dim3D() != 0)
      {
         // Check for correct sizes
         assert(k < this->mFwd.size());

         return this->mFwd.at(k).rows();
      } else
      {
         // Check for correct sizes
         assert(j < this->mFwd.size());

         return this->mFwd.at(j).rows();
      }
   }

   int TransformResolution::dimBwd(const int j, const int k) const
   {
      if(this->dim3D() != 0)
      {
         // Check for correct sizes
         assert(k < this->mBwd.size());

         return this->mBwd.at(k).rows();
      } else
      {
         // Check for correct sizes
         assert(j < this->mBwd.size());

         return this->mBwd.at(j).rows();
      }
   }

   int TransformResolution::dim2D(const int k) const
   {
      // Check for correct sizes
      assert(k < this->mIdx2D.size());

      return this->mIdx2D.at(k).size();
   }

   int TransformResolution::dim3D() const
   {
      return this->mIdx3D.size();
   }

   Datatypes::SharedScalarFieldSetup TransformResolution::spFwdSetup() const
   {
      return Datatypes::SharedScalarFieldSetup(new Datatypes::ScalarFieldSetup(this->mFwd, this->mIdx2D, this->mIdx3D));
   }

   Datatypes::SharedScalarFieldSetup TransformResolution::spBwdSetup() const
   {
      return Datatypes::SharedScalarFieldSetup(new Datatypes::ScalarFieldSetup(this->mBwd, this->mIdx2D, this->mIdx3D));
   }

   int TransformResolution::idxFwd(const int i, const int j, const int k) const
   {
      if(this->dim3D() != 0)
      {
         // Check for correct sizes
         assert(k < this->mFwd.size());
         assert(i < this->mFwd.at(k).size());

         return this->mFwd.at(k)(i);
      } else
      {
         // Check for correct sizes
         assert(j < this->mFwd.size());
         assert(i < this->mFwd.at(j).size());

         return this->mFwd.at(j)(i);
      }
   }

   int TransformResolution::idxBwd(const int i, const int j, const int k) const
   {

      if(this->dim3D() != 0)
      {
         // Check for correct sizes
         assert(k < this->mBwd.size());
         assert(i < this->mBwd.at(k).size());

         return this->mBwd.at(k)(i);
      } else
      {
         // Check for correct sizes
         assert(j < this->mBwd.size());
         assert(i < this->mBwd.at(j).size());

         return this->mBwd.at(j)(i);
      }
   }

   int TransformResolution::idx2D(const int j, const int k) const
   {
      // Check for correct sizes
      assert(k < this->mIdx2D.size());
      assert(j < this->mIdx2D.at(k).size());

      return this->mIdx2D.at(k)(j);
   }

   const ArrayI& TransformResolution::idx3D() const
   {
      return this->mIdx3D;
   }

   int TransformResolution::idx3D(const int k) const
   {
      // Check for correct sizes
      assert(k < this->mIdx3D.size());

      return this->mIdx3D(k);
   }

}
