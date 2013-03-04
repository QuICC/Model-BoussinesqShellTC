/** \file TimestepCoupling.hpp
 *  \brief Implementation of equation coupling information at the timestep level
 *
 *  \mhdBug Needs test
 */

#ifndef TIMESTEPCOUPLING_HPP
#define TIMESTEPCOUPLING_HPP

// Configuration includes
//
#include "SmartPointers/SharedPtrMacro.h"

// System includes
//

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Enums/PhysicalNames.hpp"
#include "Enums/FieldComponents.hpp"

namespace GeoMHDiSCC {

namespace Timestep {

   /**
    * @brief Implementation of equation coupling information at the timestep level
    */
   class TimestepCoupling
   {
      public:
         /// Typedef for a field ID
         typedef std::pair<PhysicalNames::Id, FieldComponents::Spectral::Id>  FieldIdType;

         /**
          * @brief Constructor
          */
        TimestepCoupling();

         /**
          * @brief Destructor
          */
         virtual ~TimestepCoupling();

         /**
          * Check if field is already present
          */
         bool isPresent(const FieldIdType& id) const;

         /**
          * @brief Get a new negative index
          */
         int newIndex();

         /**
          * @brief Get complex flag for field
          *
          * @param id   Field id
          */
         bool isComplex(const FieldIdType& id) const;

         /**
          * @brief Get index for field
          *
          * @param id   Field id
          */
         int idx(const FieldIdType& id) const;

         /**
          * @brief Update the type
          *
          * @param idx        Index to update
          * @param isComplex  Complex flag
          */
         void updateType(int idx, bool isComplex);

         /**
          * @brief Update the index
          *
          * @param oldIdx  Old index
          * @param newIdx  New index
          */
         void updateIndex(int oldIdx, int newIdx);

         /**
          * @brief Add a field
          *
          * @param id         ID of the field
          * @param isComplex  Complex flag
          * @param idx        Index
          */
         void addField(const FieldIdType& id, bool isComplex, int idx);

      protected:

      private:
         /**
          * @brief Counter for the indexes
          */
         int mCounter;

         /** 
          * Storage for the field coupling information
          */
         std::map<FieldIdType, std::pair<bool,int> >  mFields;
   };
}
}

#endif // TIMESTEPCOUPLING_HPP
