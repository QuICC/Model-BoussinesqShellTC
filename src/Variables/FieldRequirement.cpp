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
      : mIsScalar(isScalar), mNeedSpectral(needSpectral), mNeedPhysical(needPhysical), mNeedGradient(needGradient), mNeedCurl(needCurl), mPhysicalComps(3), mGradientComps(), mCurlComps(3)
   {
      // Init default physical and spectral IDs
      this->initDefaultIds();

      // Set default physical
      this->mPhysicalComps.setConstant(this->mNeedPhysical);

      // Set default curl needs
      this->mCurlComps.setConstant(this->mNeedCurl);

      // Set default physical, gradient and curl needs
      std::vector<FieldComponents::Physical::Id>::const_iterator it;
      ArrayB arr(3);
      arr.setConstant(this->mNeedGradient);
      for(it = this->mPhysicalIds.begin(); it != this->mPhysicalIds.end(); ++it)
      {
         this->mGradientComps.insert(std::make_pair(*it, arr));
      }
   }

   FieldRequirement::~FieldRequirement()
   {
   }

   void FieldRequirement::initDefaultIds()
   {
      // Set scalar spectral component
      if(this->mIsScalar)
      {
         this->mSpectralIds.push_back(FieldComponents::Spectral::SCALAR);

      // Create default spectral components
      } else
      {
         if(FieldComponents::Spectral::ONE != FieldComponents::Spectral::NOTUSED)
         {
            this->mSpectralIds.push_back(FieldComponents::Spectral::ONE);
         }
         if(FieldComponents::Spectral::TWO != FieldComponents::Spectral::NOTUSED)
         {
            this->mSpectralIds.push_back(FieldComponents::Spectral::TWO);
         }
         if(FieldComponents::Spectral::THREE != FieldComponents::Spectral::NOTUSED)
         {
            this->mSpectralIds.push_back(FieldComponents::Spectral::THREE);
         }
      }

      // Create default physical components
      if(FieldComponents::Physical::ONE != FieldComponents::Physical::NOTUSED)
      {
         this->mPhysicalIds.push_back(FieldComponents::Physical::ONE);
      }
      if(FieldComponents::Physical::TWO != FieldComponents::Physical::NOTUSED)
      {
         this->mPhysicalIds.push_back(FieldComponents::Physical::TWO);
      }
      if(FieldComponents::Physical::THREE != FieldComponents::Physical::NOTUSED)
      {
         this->mPhysicalIds.push_back(FieldComponents::Physical::THREE);
      }
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

   const ArrayB& FieldRequirement::gradientComps(const FieldComponents::Physical::Id id) const
   {
      return this->mGradientComps.find(id)->second;
   }

   const ArrayB& FieldRequirement::curlComps() const
   {
      return this->mCurlComps;
   }

   std::map<FieldComponents::Physical::Id,bool> FieldRequirement::mapPhysicalComps() const
   {
      std::map<FieldComponents::Physical::Id,bool> comps;

      for(unsigned int i = 0; i < this->mPhysicalIds.size(); i++)
      {
         comps.insert(std::make_pair(this->mPhysicalIds.at(i), this->mPhysicalComps(i)));
      }

      return comps;
   }

   std::map<FieldComponents::Physical::Id,bool> FieldRequirement::mapGradientComps(const FieldComponents::Physical::Id id) const
   {
      std::map<FieldComponents::Physical::Id,bool> comps;

      for(unsigned int i = 0; i < this->mPhysicalIds.size(); i++)
      {
         comps.insert(std::make_pair(this->mPhysicalIds.at(i), this->mGradientComps.find(id)->second(i)));
      }

      return comps;
   }

   std::map<FieldComponents::Physical::Id,bool> FieldRequirement::mapCurlComps() const
   {
      std::map<FieldComponents::Physical::Id,bool> comps;

      for(unsigned int i = 0; i < this->mPhysicalIds.size(); i++)
      {
         comps.insert(std::make_pair(this->mPhysicalIds.at(i), this->mCurlComps(i)));
      }

      return comps;
   }

   const std::vector<FieldComponents::Physical::Id>& FieldRequirement::physicalIds() const
   {
      return this->mPhysicalIds;
   }

   const std::vector<FieldComponents::Spectral::Id>& FieldRequirement::spectralIds() const
   {
      return this->mSpectralIds;
   }

   void FieldRequirement::updatePhysical(const ArrayB& comps)
   {
      this->mPhysicalComps = comps;

      this->mNeedPhysical = this->mPhysicalComps.any();
   }

   void FieldRequirement::updateGradient(const std::map<FieldComponents::Physical::Id,ArrayB>& comps)
   {
      this->mGradientComps = comps;

      // Update gradient need
      std::map<FieldComponents::Physical::Id,ArrayB>::const_iterator it;
      for(it = this->mGradientComps.begin(); it != this->mGradientComps.end(); ++it)
      {
         this->mNeedGradient = this->mNeedGradient || it->second.any();
      }
   }

   void FieldRequirement::updateCurl(const ArrayB& comps)
   {
      this->mCurlComps = comps;

      this->mNeedCurl = this->mCurlComps.any();
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
      std::vector<FieldComponents::Physical::Id>::const_iterator it;
      for(it = this->mPhysicalIds.begin(); it != this->mPhysicalIds.end(); ++it)
      {
         this->mGradientComps.find(*it)->second = this->mGradientComps.find(*it)->second + req.gradientComps(*it);
      }

      // Do OR operation of physical curl components requirement
      this->mCurlComps = this->mCurlComps + req.curlComps();
   }

}
