/** \file PeriodicOperator.hpp
 *  \brief Implementation of the multidimensional periodic operators
 */

#ifndef PERIODICOPERATOR_HPP
#define PERIODICOPERATOR_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "SpectralOperators/IOperator.hpp"

namespace GeoMHDiSCC {

namespace Spectral {

   /**
    * @brief Implementation of the multidimensional periodic operators
    */
   class PeriodicOperator
   {
      public:
         /**
          * @brief Compute the 2D Laplacian operator
          *
          * @param op   Spectral operator
          * @param k    First wave number
          * @param nBC  Number of boundary rows
          */
         static SparseMatrix laplacian2D(const IOperator& op, const MHDFloat k, const int nBC);

         /**
          * @brief Compute the 2D bi-Laplacian operator
          *
          * @param op   Spectral operator
          * @param k    First wave number
          * @param nBC  Number of boundary rows
          */
         static SparseMatrix bilaplacian2D(const IOperator& op, const MHDFloat k, const int nBC);

         /**
          * @brief Compute the 3D Laplacian operator
          *
          * @param op   Spectral operator
          * @param k    First wave number
          * @param m    Second wave number
          * @param nBC  Number of boundary rows
          */
         static SparseMatrix laplacian3D(const IOperator& op, const MHDFloat k, const MHDFloat m, const int nBC);

         /**
          * @brief Compute the 3D bi-Laplacian operator
          *
          * @param op   Spectral operator
          * @param k    First wave number
          * @param m    Second wave number
          * @param nBC  Number of boundary rows
          */
         static SparseMatrix bilaplacian3D(const IOperator& op, const MHDFloat k, const MHDFloat m, const int nBC);

         /**
          * @brief Compute the 2D quasi inverse Laplacian operator 
          *
          * @param op   Spectral operator
          * @param k    First wave number
          * @param p    Order of the quasi inverse
          */
         static SparseMatrix qLaplacian2D(const IOperator& op, const MHDFloat k, const int p);

         /**
          * @brief Compute the 2D quasi inverse bi-Laplacian operator 
          *
          * @param op   Spectral operator
          * @param k    First wave number
          * @param p    Order of the quasi inverse
          */
         static SparseMatrix qBilaplacian2D(const IOperator& op, const MHDFloat k, const int p);

         /**
          * @brief Compute the 3D quasi inverse Laplacian operator 
          *
          * @param op   Spectral operator
          * @param k    First wave number
          * @param m    Second wave number
          * @param p    Order of the quasi inverse
          */
         static SparseMatrix qLaplacian3D(const IOperator& op, const MHDFloat k, const MHDFloat m, const int p);

         /**
          * @brief Compute the 3D quasi inverse bi-Laplacian operator 
          *
          * @param op   Spectral operator
          * @param k    First wave number
          * @param m    Second wave number
          * @param p    Order of the quasi inverse
          */
         static SparseMatrix qBilaplacian3D(const IOperator& op, const MHDFloat k, const MHDFloat m, const int p);
         
      protected:

      private:
         /**
          * @brief Constructor
          */
         PeriodicOperator(const int idx, const int polyN, const ArrayI& specIdx);

         /**
          * @brief Empty Destructor
          */
         ~PeriodicOperator();
   };

}
}

#endif // PERIODICOPERATOR_HPP
