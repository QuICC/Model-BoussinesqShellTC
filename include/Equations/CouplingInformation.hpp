/** \file CouplingInformation.hpp
 *  \brief Implemenation of container the coupling information of the equations
 */

#ifndef COUPLINGINFORMATION_HPP
#define COUPLINGINFORMATION_HPP

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

namespace Equations {

   /**
    * \brief Implemenation of container the coupling information of the equations
    */
   class CouplingInformation
   {
      public:
         /// Typedef to simplify notation for the field coupling data
         typedef std::multimap<FieldComponents::Spectral::Id, std::pair<PhysicalNames::Id, FieldComponents::Spectral::Id> > FieldCouplingType;

         /// Typedef for an iterator for the field coupling data
         typedef FieldCouplingType::const_iterator  field_iterator;

         /// Typedef for a range iterator for the field coupling data
         typedef std::pair<field_iterator,field_iterator>  field_iterator_range;

         /**
          * @brief Simple constructor
          */
         CouplingInformation();

         /**
          * @brief Simple empty destructor
          */
         virtual ~CouplingInformation();

         /**
          * @brief Set a field coupling
          *
          * @param comp    Field component
          * @param field   External coupling information
          */
         void addField(const FieldComponents::Spectral::Id& comp, const std::pair<PhysicalNames::Id, FieldComponents::Spectral::Id>& field);

         /**
          * @brief Set an internal coupling (matrix dimensions)
          *
          * @param comp Component ID
          * @param nMat Number of matrices
          * @param dim  Number of columns on RHS
          */
         void addInternal(const FieldComponents::Spectral::Id& comp, const int nMat, const ArrayI& dim);

         /**
          * @brief Get internal coupling information (matrix dimensions)
          */
         const std::pair<int,ArrayI>&  internal(const FieldComponents::Spectral::Id& comp) const;

         /**
          * @brief Equation has external coupling (i.e. with other fields)?
          */
         bool hasFieldCoupling() const;

         /**
          * @brief Get the number of field couplings
          */
         int nFields(const FieldComponents::Spectral::Id& comp) const;

         /**
          * @brief Get iterator to field couplings
          */
         field_iterator_range fieldRange(const FieldComponents::Spectral::Id& comp) const;

      protected:

      private:
         /**
          * @brief Other field coupling flag
          */
         bool mHasField;

         /**
          * @brief Internal coupling flag
          */
         bool mHasInternal;

         /**
          * @brief Storage for the field coupling information
          */
         std::multimap<FieldComponents::Spectral::Id, std::pair<PhysicalNames::Id, FieldComponents::Spectral::Id> >   mField;

         /**
          * @brief Storage for the internal coupling information
          */
         std::map<FieldComponents::Spectral::Id, std::pair<int, ArrayI> >   mInternal;
   };

   inline bool CouplingInformation::hasFieldCoupling() const
   {
      return this->mHasField;
   }

   inline const std::pair<int,ArrayI>& CouplingInformation::internal(const FieldComponents::Spectral::Id& comp) const
   {
      return this->mInternal.find(comp)->second;
   }

   inline int CouplingInformation::nFields(const FieldComponents::Spectral::Id& comp) const
   {
      return this->mField.count(comp);
   }

   inline CouplingInformation::field_iterator_range CouplingInformation::fieldRange(const FieldComponents::Spectral::Id& comp) const
   {
      return this->mField.equal_range(comp);
   }
}
}

#endif // COUPLINGINFORMATION_HPP
