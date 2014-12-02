/** 
 * @file FieldRequirement.cpp
 * @brief Source of the variable requirements
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "Variables/FieldRequirement.hpp"

// Project includes
//

namespace GeoMHDiSCC {

   FieldRequirement::FieldRequirement(const bool isScalar, const bool needSpectral, const bool needPhysical, const bool needGradient, const bool needCurl)
      : mIsScalar(isScalar), mNeedSpectral(needSpectral), mNeedPhysical(needPhysical), mNeedGradient(needGradient), mNeedCurl(needCurl), mPhysicalComps(3), mGradientComps(3), mCurlComps(3)
   {
      this->mPhysicalComps.setConstant(this->mNeedPhysical);
      this->mGradientComps.setConstant(this->mNeedGradient);
      this->mCurlComps.setConstant(this->mNeedCurl);
   }

   FieldRequirement::FieldRequirement(const bool isScalar, const bool needSpectral, const ArrayB& physicalComps, const ArrayB& gradientComps, const ArrayB& curlComps)
      : mIsScalar(isScalar), mNeedSpectral(needSpectral), mNeedPhysical(physicalComps.any()), mNeedGradient(gradientComps.any()), mNeedCurl(curlComps.any()), mPhysicalComps(physicalComps), mGradientComps(gradientComps), mCurlComps(curlComps)
   {
   }

   FieldRequirement::~FieldRequirement()
   {
   }

   bool FieldRequirement::isScalar() const
   {
      return this->mIsScalar;
   }

   bool FieldRequirement::needSpectral() const
   {
      return this->mNeedSpectral;
   }

   bool FieldRequirement::needPhysical() const
   {
      return this->mNeedPhysical;
   }

   bool FieldRequirement::needPhysicalGradient() const
   {
      return this->mNeedGradient;
   }

   bool FieldRequirement::needPhysicalCurl() const
   {
      return this->mNeedCurl;
   }

   const ArrayB& FieldRequirement::physicalComps() const
   {
      return this->mPhysicalComps;
   }

   const ArrayB& FieldRequirement::gradientComps() const
   {
      return this->mGradientComps;
   }

   const ArrayB& FieldRequirement::curlComps() const
   {
      return this->mCurlComps;
   }

   void FieldRequirement::merge(const FieldRequirement& req)
   {
      // Assert for same type
      assert(this->mIsScalar == req.isScalar());

      // Do OR operation on spectral requirement
      this->mNeedSpectral = this->mNeedSpectral || req.needSpectral(); 

      // Do OR operation on physical requirement
      this->mNeedPhysical = this->mNeedPhysical || req.needPhysical(); 

      // Do OR operation on physical gradient requirement
      this->mNeedGradient = this->mNeedGradient || req.needPhysicalGradient(); 

      // Do OR operation on physical curl requirement
      this->mNeedCurl = this->mNeedCurl || req.needPhysicalCurl(); 

      // Do OR operation of physical components requirement
      this->mPhysicalComps = this->mPhysicalComps + req.physicalComps();

      // Do OR operation of physical gradient components requirement
      this->mGradientComps = this->mGradientComps + req.gradientComps();

      // Do OR operation of physical curl components requirement
      this->mCurlComps = this->mCurlComps + req.curlComps();
   }

}
