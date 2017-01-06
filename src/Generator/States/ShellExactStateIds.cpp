/** 
 * @file ShellExactStateIds.cpp
 * @brief Source of the implementation of the tools to generate exact spherical shell states
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Configuration includes
//

// System includes
//
#include <tr1/cmath>

// External includes
//

// Class include
//
#include "Generator/States/ShellExactStateIds.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"

namespace QuICC {

namespace Equations {

	const MHDFloat ShellExactStateIds::PCOS = 99999;

   const MHDFloat ShellExactStateIds::PSIN = -99999;

   Array ShellExactStateIds::sph_harmonic(const MHDComplex amplitude, const int l, const int m, const MHDFloat theta, const Array& phi)
   {
      MHDFloat re = amplitude.real();
      MHDFloat im = amplitude.imag();

      Array sph_harm = 2.0*re*(static_cast<MHDFloat>(m)*phi).array().cos() - 2.0*im*(static_cast<MHDFloat>(m)*phi).array().sin();
      #if defined QUICC_SHNORM_SCHMIDT
         // Schmidt quasi-normalized spherical harmonic Y_l^m
         MHDFloat leg = std::tr1::sph_legendre(l, m, theta)*std::sqrt(4.0*Math::PI/static_cast<MHDFloat>(2*l + 1));
      #elif defined QUICC_SHNORM_UNITY
         // Normalized spherical harmonic Y_l^m
         MHDFloat leg = std::tr1::sph_legendre(l, m, theta);
      #endif //defined QUICC_SHNORM_SCHMIDT
      
      sph_harm *= leg;

      return sph_harm;
   }

}
}
