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
#include "Enums/FieldComponents.hpp"
#include "Enums/PhysicalNames.hpp"

namespace GeoMHDiSCC {

namespace Timestep {

   /**
    * \brief Implementation of the base for an equation timestepper
    */
   class EquationTimestepperBase
   {
      public:
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
         void addInformation(const std::pair<PhysicalNames::Id, FieldComponents::Spectral::Id>& id, const ArrayI& startRow);

         /**
          * @brief Get start row 
          */
         int startRow(const std::pair<PhysicalNames::Id, FieldComponents::Spectral::Id>& id, const int i) const;
         
      protected:
         /**
          * @brief Starting index
          */
         int mZeroIdx;

         /**
          * @brief Storage for the storage information
          */
         std::map<std::pair<PhysicalNames::Id, FieldComponents::Spectral::Id>, ArrayI> mInformation;

      private:
   };

   /// Typedef for a shared pointer of a EquationTimestepperBase
   typedef SharedPtrMacro<EquationTimestepperBase>  SharedEquationTimestepperBase;
}
}

#endif // EQUATIONTIMESTEPPERBASE_HPP
