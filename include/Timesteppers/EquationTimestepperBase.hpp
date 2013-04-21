/** \file EquationTimestepperBase.hpp
 *  \brief Implementation of the base for an equation timestepper
 */

#ifndef EQUATIONTIMESTEPPERBASE_HPP
#define EQUATIONTIMESTEPPERBASE_HPP

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
#include "Enums/FieldIds.hpp"

namespace GeoMHDiSCC {

namespace Timestep {

   /**
    * \brief Implementation of the base for an equation timestepper
    */
   class EquationTimestepperBase
   {
      public:
         /// Typedef to simplify notation for the field data
         typedef std::vector<SpectralFieldId> FieldIdVector;

         /// Typedef for an iterator for the field data
         typedef FieldIdVector::const_iterator  field_iterator;

         /// Typedef for a range iterator for the field coupling data
         typedef std::pair<field_iterator,field_iterator>  field_iterator_range;

         /**
          * @brief Constructor
          *
          * @param start   Starting index (for example without m=0)
          */
         EquationTimestepperBase(const int start);

         /**
          * @brief Destructor
          */
         virtual ~EquationTimestepperBase();

         /**
          * @brief Add storage information 
          */
         void addInformation(const SpectralFieldId& id, const ArrayI& startRow);

         /**
          * @brief Get start row 
          */
         int startRow(const SpectralFieldId& id, const int i) const;

         /**
          * @brief Range of stored fields
          */
         field_iterator_range fieldRange() const;
         
      protected:
         /**
          * @brief Starting index
          */
         int mZeroIdx;

         /**
          * @brief Storage for the field Ids
          */
         std::vector<SpectralFieldId>   mFieldIds;

         /**
          * @brief Storage for the storage information
          */
         std::map<SpectralFieldId, ArrayI> mInformation;

      private:
   };

   /// Typedef for a shared pointer of a EquationTimestepperBase
   typedef SharedPtrMacro<EquationTimestepperBase>  SharedEquationTimestepperBase;
}
}

#endif // EQUATIONTIMESTEPPERBASE_HPP
