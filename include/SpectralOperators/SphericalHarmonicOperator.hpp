/** \file SphericalHarmonicOperator.hpp
 *  \brief Implementation of the spherical harmonic operators
 */

#ifndef SPHERICALHARMONICOPERATOR_HPP
#define SPHERICALHARMONICOPERATOR_HPP

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
    * @brief Implementation of the spherical harmonic operators
    */
   class SphericalHarmonicOperator
   {
      public:
         /**
          * @brief Compute the spherical surface Laplacian operator
          *
          * @param l    Harmonic degree
          * @param m    Harmonic order
          */
         static MHDFloat surfaceLaplacian(const MHDFloat l, const MHDFloat m);

         /**
          * @brief Compute the spherical surface bi-Laplacian
          *
          * @param l    Harmonic degree
          * @param m    Harmonic order
          */
         static MHDFloat surfaceBilaplacian(const MHDFloat l, const MHDFloat m);

         /**
          * @brief Compute the spherical Laplacian operator
          *
          * @param op   Radial spectral operator
          * @param l    Harmonic degree
          * @param m    Harmonic order
          * @param nBC  Number of boundary rows
          */
         static SparseMatrix laplacian(const IOperator& op, const MHDFloat l, const MHDFloat m, const int nBC);

         /**
          * @brief Compute the spherical bi-Laplacian operator
          *
          * @param op   Radial spectral operator
          * @param l    Harmonic degree
          * @param m    Harmonic order
          * @param nBC  Number of boundary rows
          */
         static SparseMatrix bilaplacian(const IOperator& op, const MHDFloat k, const MHDFloat m, const int nBC);

         /**
          * @brief Compute the quasi inverse spherical Laplacian operator 
          *
          * @param op   Radial spectral operator
          * @param l    Harmonic degree
          * @param m    Harmonic order
          * @param p    Order of the quasi inverse
          */
         static SparseMatrix qLaplacian(const IOperator& op, const MHDFloat l, const MHDFloat m, const int p);

         /**
          * @brief Compute the quasi inverse spherical bi-Laplacian operator 
          *
          * @param op   Radial spectral operator
          * @param l    Harmonic degree
          * @param m    Harmonic order
          * @param p    Order of the quasi inverse
          */
         static SparseMatrix qBilaplacian(const IOperator& op, const MHDFloat l, const MHDFloat m, const int p);
         
      protected:

      private:
         /**
          * @brief Constructor
          */
         SphericalHarmonicOperator();

         /**
          * @brief Empty Destructor
          */
         ~SphericalHarmonicOperator();
   };

}
}

#endif // SPHERICALHARMONICOPERATOR_HPP
