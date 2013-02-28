/** \file EquationTimestepperBase.hpp
 *  \brief Implementation of the base for an equation timestepper
 *
 *  \mhdBug Needs test
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
#include "Enums/FieldComponents.hpp"
#include "Enums/PhysicalNames.hpp"

namespace GeoMHDiSCC {

   /**
    * \brief Implementation of the base for an equation timestepper
    */
   class EquationTimestepperBase
   {
      public:
         /**
          * @brief Constructor
          *
          * @param nField Number of fields
          */
         EquationTimestepperBase(const int nField);

         /**
          * @brief Destructor
          */
         virtual ~EquationTimestepperBase();

         /**
          * @brief Get index of current field
          */
         int current() const;

         /**
          * @brief Move counter to next field (loops around)
          */
         void next();

         /**
          * @brief Add storage information 
          */
         void addInformation(const std::pair<PhysicalNames::Id, FieldComponents::Spectral::Id>& id, const ArrayI& startRow);

         /**
          * @brief Get start row 
          */
         int startRow(const std::pair<PhysicalNames::Id, FieldComponents::Spectral::Id>& id, const int i) const;
         
      protected:
         /**
          * @brief Storage for the number of fields
          */
         int mNField;

         /**
          * @brief Storage for the current field index
          */
         int mCurrent;

         /**
          * @brief Storage for the storage information
          */
         std::map<std::pair<PhysicalNames::Id, FieldComponents::Spectral::Id>, ArrayI> mInformation;

      private:
   };

   inline int EquationTimestepperBase::current() const
   {
      return this->mCurrent;
   }

   inline void EquationTimestepperBase::next()
   {
      this->mCurrent = (this->mCurrent+1) % this->mNField;
   }

   /// Typedef for a shared pointer of a EquationTimestepperBase
   typedef SharedPtrMacro<EquationTimestepperBase>  SharedEquationTimestepperBase;
}

#endif // EQUATIONTIMESTEPPERBASE_HPP
