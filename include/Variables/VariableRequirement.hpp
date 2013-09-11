/** 
 * @file VariableRequirement.hpp
 * @brief Implementation of a class to store requirements for the variables
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 * @version 0.9.0
 * @date 2013-09-11
 */

#ifndef VARIABLEREQUIREMENT_HPP
#define VARIABLEREQUIREMENT_HPP

// System includes
//
#include <map>

// External includes
//

// Project includes
//
#include "Enums/FieldIds.hpp"
#include "Variables/FieldRequirement.hpp"

namespace GeoMHDiSCC {

   /**
    * @brief Implementation of a class to store requirements for the variables
    */
   class VariableRequirement
   {
      public:
         /// Typedef for the const iterator
         typedef  std::map<PhysicalNames::Id,FieldRequirement>::const_iterator   const_iterator;

         /**
          * @brief Constructor
          */
         VariableRequirement();

         /**
          * @brief Destructor
          */
         ~VariableRequirement();

         /**
          * @brief Get field requirements
          */
         const FieldRequirement& field(const PhysicalNames::Id id) const;

         /**
          * @brief Add field requirement
          */
         void addField(const PhysicalNames::Id id, const FieldRequirement& req);

         /**
          * @brief Const iterator to access all requirements
          */
         const_iterator begin() const;

         /**
          * @brief Const iterator to access all requirements
          */
         const_iterator end() const;

         /**
          * @brief Merge the information of two requirements
          *
          * The operation is basically a OR on the need? calls
          */
         void merge(const VariableRequirement& req);
         
      protected:

      private:
         /**
          * @brief Not required field (all answers are false)
          */
         FieldRequirement  mNoField;

         /**
          * @brief Storage for all the information
          */
         std::map<PhysicalNames::Id, FieldRequirement> mInfo;
   };

}

#endif // VARIABLEREQUIREMENT_HPP
