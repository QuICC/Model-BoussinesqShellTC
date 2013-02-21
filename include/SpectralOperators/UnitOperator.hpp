/** \file UnitOperator.hpp
 *  \brief Implementation of the unit spectral operator with quasi inverses
 */

#ifndef UNITOPERATOR_HPP
#define UNITOPERATOR_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "SpectralOperators/IOperator.hpp"
#include "Simulation/Enums/BoundaryConditions.hpp"

namespace GeoMHDiSCC {

namespace Spectral {

   /**
    * \brief Implementation of the unit spectral operator with quasi inverses
    */
   class UnitOperator: IOperator
   {
      public:
         /**
          * @brief Constructor
          *
          * @param idx     Index of assigned parameter
          * @param polyN   Size of the polynomial basis
          * @param specIdx Spectral indexes
          */
         UnitOperator(const int idx, const int polyN, const ArrayI& specIdx);

         /**
          * @brief Empty Destructor
          */
         ~UnitOperator();

         /**
          * @brief (Possibly) change dimension and spectral indexes. Returns true if new operators a created
          *
          * This routine is to be used when looping over all indexes to extract the minimal set of required matrices.
          *
          * @param polyN   New dimension
          * @param specIdx New set of spectral indexes
          */
         bool loopNext(const int polyN, const ArrayI& specIdx);

         /**
          * @brief Get the identity matrix
          *
          * @param p Order of the quasi identity
          */
         SparseMatrix id(const int p = 0);

         /**
          * @brief Get the derivative of order p
          *
          * @param nBC  Number of boundary
          * @param p Order of the derivative
          */
         SparseMatrix diff(const int nBC, const int p);

         /**
          * @brief Get the quasi inverse derivative of order p
          *
          * @param p    Order of the quasi inverse
          * @param q    Order of the derivative
          */
         SparseMatrix qDiff(const int p, const int q);

         /**
          * @brief Convert tau lines to complete sparse matrix
          *
          * @param bcId List of BC to implement
          */
         DecoupledZSparse tau(const std::map<BoundaryConditions::Id,BoundaryConditions::Position>& bcId);
         
      protected:

      private:
   };

}
}

#endif // UNITOPERATOR_HPP
