/** 
 * @file TauChebyshev.hpp
 * @brief Implementation of the Tau conditions for chebyshev based operators
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef TAUCHEBYSHEV_HPP
#define TAUCHEBYSHEV_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"

namespace GeoMHDiSCC {

namespace Spectral {

   /**
    * @brief Implementation of the Tau conditions for chebyshev based operators
    */
   class TauChebyshev
   {
      public:
         /**
          * @brief Constructor
          *
          * @param nN   Size of the Tau basis
          * @param bcs  Vector of boundary conditions
          * @param nEq  Number of equations (only independent of BC for coupled systems)
          */
         TauChebyshev(const int nN, const Boundary::BCvector& bcs, const int nEq);

         /**
          * @brief Empty Destructor
          */
         ~TauChebyshev();

         /**
          * @brief Constrain matrix with boundary conditions
          *
          * @param mat  Matrix operator to constrain
          */
         template <typename TData> const TData& constrain(const TData& mat);

         /**
          * @brief Extend Galerkin basis coefficients to tau expansion 
          *
          * @param spec Matrix of spectral coefficients
          */
         template <typename TData> const TData& extend(const TData& spec);

         /**
          * @brief Restrict tau coefficients to galerkin size 
          *
          * @param spec Matrix of spectral coefficients
          */
         template <typename TData> const TData& restrict(const TData& spec);
         
      protected:

      private:
   };

   template <typename TData> inline const TData& TauChebyshev::extend(const TData& spec)
   {
      return spec;
   }

   template <typename TData> inline const TData& TauChebyshev::restrict(const TData& spec)
   {
      return spec;
   }

   template <typename TData> const TData& TauChebyshev::constrain(const TData& mat)
   {
   }

}
}

#endif // TAUCHEBYSHEV_HPP
