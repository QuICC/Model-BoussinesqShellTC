/** 
 * @file FieldRequirement.hpp
 * @brief Implementation of a class to store requirements for a field
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef FIELDREQUIREMENT_HPP
#define FIELDREQUIREMENT_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"

namespace GeoMHDiSCC {

   /**
    * @brief Implementation of a class to store requirements for a field
    */
   class FieldRequirement 
   {
      public:
         /**
          * @brief General constructor without detailed field components
          */
         FieldRequirement(const bool isScalar, const bool needSpectral, const bool needPhysical, const bool needGradient, const bool needCurl = false);

         /**
          * @brief Constructor requiring detailed field components requirements
          */
         FieldRequirement(const bool isScalar, const bool needSpectral, const ArrayB& physicalComps, const ArrayB& gradientComps, const ArrayB& curlComps);

         /**
          * @brief Destructor
          */
         ~FieldRequirement();

         /**
          * @brief Check if it is scalar field
          */
         bool isScalar() const;

         /**
          * @brief Get spectral requirement
          */
         bool needSpectral() const;

         /**
          * @brief Get physical field requirement
          */
         bool needPhysical() const;

         /**
          * @brief Get physical gradient requirement
          */
         bool needPhysicalGradient() const;

         /**
          * @brief Get physical curl requirement
          */
         bool needPhysicalCurl() const;

         /**
          * @brief Get the physical field components requirements
          */
         const ArrayB& physicalComps() const;

         /**
          * @brief Get the physical gradient components requirements
          */
         const ArrayB& gradientComps() const;

         /**
          * @brief Get the physical curl components requirements
          */
         const ArrayB& curlComps() const;

         /**
          * @brief Merge information from other requirements
          *
          * @param req Requirements to merge
          */
         void merge(const FieldRequirement& req);
         
      protected:

      private:
         /**
          * @brief Is field a scalar field?
          */
         bool  mIsScalar;

         /**
          * @brief Storage for spectral storage requirements
          */
         bool  mNeedSpectral;

         /**
          * @brief Storage for physical storage requirements
          */
         bool  mNeedPhysical;

         /**
          * @brief Storage for physical gradient storage requirements
          */
         bool  mNeedGradient;

         /**
          * @brief Storage for physical curl storage requirements
          */
         bool  mNeedCurl;

         /**
          * @brief Detailed requirements for physical field components
          */
         ArrayB mPhysicalComps;

         /**
          * @brief Detailed requirements for gradient field components
          */
         ArrayB mGradientComps;

         /**
          * @brief Detailed requirements for curl field components
          */
         ArrayB mCurlComps;
   };

}

#endif // FIELDREQUIREMENT_HPP
