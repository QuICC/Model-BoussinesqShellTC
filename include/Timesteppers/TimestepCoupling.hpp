/** \file TimestepCoupling.hpp
 *  \brief Implementation of equation coupling information at the timestep level
 */

#ifndef TIMESTEPCOUPLING_HPP
#define TIMESTEPCOUPLING_HPP

// Configuration includes
//
#include "SmartPointers/SharedPtrMacro.h"

// System includes
//
#include <tr1/tuple>

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
         void updateType(const int idx, const bool isComplex);

         /**
          * @brief Check starting index
          *
          * @param idx     Index to update
          * @param start   Starting index
          */
         void checkStart(const int idx, const int start);

         /**
          * @brief Update the index
          *
          * @param oldIdx  Old index
          * @param newIdx  New index
          */
         void updateIndex(const int oldIdx, const int newIdx);

         /**
          * @brief Add a field
          *
          * @param id         ID of the field
          * @param isComplex  Complex flag
          * @param idx        Index
          * @param start      Starting index
          */
         void addField(const FieldIdType& id, const bool isComplex, const int idx, const int start);

      protected:

      private:
         /**
          * @brief Counter for the indexes
          */
         int mCounter;

         /** 
          * Storage for the field coupling information
          */
         std::map<FieldIdType, std::tr1::tuple<bool,int,int> >  mFields;
   };
}
}

#endif // TIMESTEPCOUPLING_HPP
