/** \file FieldRequirement.hpp
 *  \brief Implementation of a class to store requirements for a field
 *
 *  \mhdBug Needs test
 */

#ifndef FIELDREQUIREMENT_HPP
#define FIELDREQUIREMENT_HPP

// System includes
//

// External includes
//

// Project includes
//

namespace GeoMHDiSCC {

   /**
    * @brief Implementation of a class to store requirements for a field
    */
   class FieldRequirement 
   {
      public:
         /**
          * @brief Constructor
          */
         FieldRequirement(const bool isScalar, const bool needSpectral, const bool needPhysical, const bool needDiff);

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
          * @brief Get spectral requirement
          */
         bool needPhysical() const;

         /**
          * @brief Get spectral requirement
          */
         bool needPhysicalDiff() const;

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
          * @brief Storage for physical differential storage requirements
          */
         bool  mNeedDiff;
   };

}

#endif // FIELDREQUIREMENT_HPP
