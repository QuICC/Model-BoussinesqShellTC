/** \file ImExRK3.cpp
 *  \brief Implementation of an implicit/explicit Runge-Kutta scheme of order ~3
 */

// System includes
//

// External includes
//

// Class include
//
#include "Timesteppers/ImExRK3.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Timestep {

   // Scheme requires 3 substeps
   const int ImExRK3::STEPS = 3;

   // Set the alpha parameters
   const Eigen::Array<MHDFloat,3,1> ImExRK3::mAlpha = Eigen::Array<MHDFloat,3,1>(29./96., -3./40., 1./6.);

   // Set the beta parameters
   const Eigen::Array<MHDFloat,3,1> ImExRK3::mBeta = Eigen::Array<MHDFloat,3,1>(37./160., 5./24., 1./6.);

   // Set the gamma parameters
   const Eigen::Array<MHDFloat,3,1> ImExRK3::mGamma = Eigen::Array<MHDFloat,3,1>(8./15., 5./12., 3./4.);

   // Set the zeta parameters
   const Eigen::Array<MHDFloat,3,1> ImExRK3::mZeta = Eigen::Array<MHDFloat,3,1>(0., -17./60., -5./12.);

   MHDFloat ImExRK3::lhsT(const int step)
   {
      return 1./ImExRK3::mBeta(step);
   }

   MHDFloat ImExRK3::lhsL(const int step)
   {
      return 1.0;
   }

   MHDFloat ImExRK3::rhsT(const int step)
   {
      return 1./ImExRK3::mBeta(step);
   }

   MHDFloat ImExRK3::rhsL(const int step)
   {
      return ImExRK3::mAlpha(step)/ImExRK3::mBeta(step);
   }

   MHDFloat ImExRK3::rhsN(const int step)
   {
      return ImExRK3::mGamma(step)/ImExRK3::mBeta(step);
   }

   MHDFloat ImExRK3::rhsNN(const int step)
   {
      return ImExRK3::mZeta(step)/ImExRK3::mBeta(step);
   }

}
}
