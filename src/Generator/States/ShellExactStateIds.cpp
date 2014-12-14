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

namespace GeoMHDiSCC {

namespace Equations {

   Array ShellExactStateIds::sph_harmonic(const MHDComplex amplitude, const int l, const int m, const MHDFloat theta, const Array& phi)
   {
      MHDFloat re = amplitude.real();
      MHDFloat im = amplitude.imag();

      Array sph_harm = re*(static_cast<MHDFloat>(m)*phi).array().cos() + im*(static_cast<MHDFloat>(m)*phi).array().sin();
      MHDFloat leg = std::tr1::sph_legendre(l, m, theta);
      sph_harm *= leg;

      return sph_harm;
   }

}
}
