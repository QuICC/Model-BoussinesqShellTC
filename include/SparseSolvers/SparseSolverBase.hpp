/** \file SparseSolverBase.hpp
 *  \brief Implementation of the base for the solver structures
 */

#ifndef SPARSESOLVERBASE_HPP
#define SPARSESOLVERBASE_HPP

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

namespace Solver {

   /**
    * \brief Implementation of the base for the solver structures
    */
   class SparseSolverBase
   {
      public:
         /// Typedef to simplify notation for the field data
         typedef std::vector<SpectralFieldId> FieldIdVector;

         /// Typedef for an iterator for the field data
         typedef FieldIdVector::const_iterator  FieldId_iterator;

         /// Typedef for a range iterator for the field coupling data
         typedef std::pair<FieldId_iterator,FieldId_iterator>  FieldId_range;

         /**
          * @brief Constructor
          *
          * @param start   Starting index (for example without m=0)
          */
         SparseSolverBase(const int start);

         /**
          * @brief Destructor
          */
         virtual ~SparseSolverBase();

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
         FieldId_range fieldRange() const;
         
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

   /// Typedef for a shared pointer of a SparseSolverBase
   typedef SharedPtrMacro<SparseSolverBase>  SharedSparseSolverBase;
}
}

#endif // SPARSESOLVERBASE_HPP
