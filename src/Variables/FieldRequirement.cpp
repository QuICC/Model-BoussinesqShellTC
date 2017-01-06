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
#include "Exceptions/Exception.hpp"
#include "Variables/FieldRequirement.hpp"

// Project includes
//

namespace QuICC {

   FieldRequirement::FieldRequirement(const bool isScalar, const bool needSpectral, const bool needPhysical, const bool needGradient, const bool needCurl, const bool needGradient2)
      : mIsScalar(isScalar), mNeedSpectral(needSpectral), mNeedPhysical(needPhysical), mNeedGradient(needGradient), mNeedCurl(needCurl), mNeedGradient2(needGradient2), mPhysicalComps(3), mGradientComps(), mCurlComps(3), mGradient2Comps()
   {
      // Stop if curl is requested for scalar
      if(isScalar && needCurl)
      {
         throw Exception("Tried to setup curl calculation for a scalar!");
      }

      // Init default physical and spectral IDs
      this->initDefaultIds();

      // Set default physical
      this->mPhysicalComps.setConstant(this->mNeedPhysical);

      // Set default curl needs
      this->mCurlComps.setConstant(this->mNeedCurl);

      // Set default physical, gradient and curl needs
      std::vector<FieldComponents::Spectral::Id>::const_iterator it;
      ArrayB arr(3);
      arr.setConstant(this->mNeedGradient);
      for(it = this->mSpectralIds.begin(); it != this->mSpectralIds.end(); ++it)
      {
         this->mGradientComps.insert(std::make_pair(*it, arr));
      }

      // Set default 2nd order gradient needs
      MatrixB mat = MatrixB::Zero(3,3);
      mat.triangularView<Eigen::Upper>().setConstant(this->mNeedGradient2);
      for(it = this->mSpectralIds.begin(); it != this->mSpectralIds.end(); ++it)
      {
         this->mGradient2Comps.insert(std::make_pair(*it, mat));
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

   bool FieldRequirement::needPhysicalGradient2() const
   {
      return this->mNeedGradient2;
   }

   const ArrayB& FieldRequirement::physicalComps() const
   {
      return this->mPhysicalComps;
   }

   const ArrayB& FieldRequirement::gradientComps(const FieldComponents::Spectral::Id id) const
   {
      assert(this->mGradientComps.count(id));

      return this->mGradientComps.find(id)->second;
   }

   const ArrayB& FieldRequirement::curlComps() const
   {
      return this->mCurlComps;
   }

   const MatrixB& FieldRequirement::gradient2Comps(const FieldComponents::Spectral::Id id) const
   {
      assert(this->mGradient2Comps.count(id));

      return this->mGradient2Comps.find(id)->second;
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

   std::map<FieldComponents::Physical::Id,bool> FieldRequirement::mapGradientComps(const FieldComponents::Spectral::Id id) const
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

   std::map<std::pair<FieldComponents::Physical::Id,FieldComponents::Physical::Id>,bool> FieldRequirement::mapGradient2Comps(const FieldComponents::Spectral::Id id) const
   {
      std::map<std::pair<FieldComponents::Physical::Id,FieldComponents::Physical::Id>,bool> comps;

      for(unsigned int i = 0; i < this->mPhysicalIds.size(); i++)
      {
         for(unsigned int j = 0; j < this->mPhysicalIds.size(); j++)
         {
            comps.insert(std::make_pair(std::make_pair(this->mPhysicalIds.at(i),this->mPhysicalIds.at(j)), this->mGradient2Comps.find(id)->second(i,j)));
         }
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

   void FieldRequirement::updateGradient(const std::map<FieldComponents::Spectral::Id,ArrayB>& comps)
   {
      this->mGradientComps = comps;

      // Update gradient need
      std::map<FieldComponents::Spectral::Id,ArrayB>::const_iterator it;
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

   void FieldRequirement::updateGradient2(const std::map<FieldComponents::Spectral::Id,MatrixB>& comps)
   {
      this->mGradient2Comps = comps;

      // Update 2nd order gradient need
      std::map<FieldComponents::Spectral::Id,MatrixB>::iterator it;
      for(it = this->mGradient2Comps.begin(); it != this->mGradient2Comps.end(); ++it)
      {
         // Clear strictly lower triangular part
         it->second.triangularView<Eigen::StrictlyLower>().setZero();

         this->mNeedGradient2 = this->mNeedGradient2 || it->second.any();
      }
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

      // Do OR operation on physical 2nd order gradient requirement
      this->mNeedGradient2 = this->mNeedGradient2 || req.needPhysicalGradient2(); 

      // Do OR operation of physical components requirement
      this->mPhysicalComps = this->mPhysicalComps + req.physicalComps();

      // Do OR operation of physical gradient components requirement
      if(req.needPhysicalGradient())
      {
         std::vector<FieldComponents::Spectral::Id>::const_iterator it;
         for(it = this->mSpectralIds.begin(); it != this->mSpectralIds.end(); ++it)
         {
            this->mGradientComps.find(*it)->second = this->mGradientComps.find(*it)->second + req.gradientComps(*it);
         }
      }

      // Do OR operation of physical curl components requirement
      if(req.needPhysicalCurl())
      {
         this->mCurlComps = this->mCurlComps + req.curlComps();
      }

      // Do OR operation of physical 2nd order gradient components requirement
      if(req.needPhysicalGradient2())
      {
         std::vector<FieldComponents::Spectral::Id>::const_iterator it;
         for(it = this->mSpectralIds.begin(); it != this->mSpectralIds.end(); ++it)
         {
            this->mGradient2Comps.find(*it)->second = this->mGradient2Comps.find(*it)->second + req.gradient2Comps(*it);
         }
      }
   }

}
