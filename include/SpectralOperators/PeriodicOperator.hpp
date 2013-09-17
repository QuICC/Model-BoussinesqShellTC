/** 
 * @file PeriodicOperator.hpp
 * @brief Implementation of the multidimensional periodic operators
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
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
          * @brief Compute the 2D Laplacian operator for a periodic 2D box
          *
          * @param k    First wave number
          * @param m    Second wave number
          */
         static MHDFloat laplacian2D(const MHDFloat k, const MHDFloat m);

         /**
          * @brief Compute the 2D bi-Laplacian operator for a periodic 2D box
          *
          * @param k    First wave number
          * @param m    Second wave number
          */
         static MHDFloat bilaplacian2D(const MHDFloat k, const MHDFloat m);

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
          * @param opA  First spectral operator
          * @param opB  Second spectral operator
          * @param k    First wave number
          * @param nBCA  Number of boundary rows for first operator
          * @param nBCB  Number of boundary rows for second operator
          */
         static SparseMatrix laplacian3D(const IOperator& opA, const IOperator& opB, const MHDFloat k, const int nBCA, const int nBCB);

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
          * @param opA  First spectral operator
          * @param opB  Second spectral operator
          * @param k    First wave number
          * @param pA   Order of the quasi inverse for first operator
          * @param pB   Order of the quasi inverse for second operator
          */
         static SparseMatrix qLaplacian3D(const IOperator& opA, const IOperator& opB, const MHDFloat k, const int pA, const int pB);

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
          * @param opA  First spectral operator
          * @param opB  Second spectral operator
          * @param k    First wave number
          * @param pA   Order of the quasi inverse for first operator
          * @param pB   Order of the quasi inverse for second operator
          */
         static SparseMatrix qBilaplacian3D(const IOperator& opA, const IOperator& opB, const MHDFloat k, const int pA, const int pB);

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
         PeriodicOperator();

         /**
          * @brief Empty Destructor
          */
         ~PeriodicOperator();
   };

}
}

#endif // PERIODICOPERATOR_HPP
