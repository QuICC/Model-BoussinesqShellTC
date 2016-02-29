/** 
 * @file ImExRK3.cpp
 * @brief Implementation of an implicit/explicit Runge-Kutta scheme of order ~3
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
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

   // Scheme order
   const int ImExRK3::ORDER = 2;

   // Scheme requires field at t_n
   const int ImExRK3::FIELD_MEMORY = 0;

   // Scheme requires nonlinear term at t(n-1)
   const int ImExRK3::NONLINEAR_MEMORY = 1;

   // Name of the scheme
   const std::string ImExRK3::NAME = "ImExRK3";

   // Set the alpha parameters
   const Eigen::Array<MHDFloat,3,1> ImExRK3::mAlpha = Eigen::Array<MHDFloat,3,1>(29./96., -3./40., 1./6.);

   // Set the beta parameters
   const Eigen::Array<MHDFloat,3,1> ImExRK3::mBeta = Eigen::Array<MHDFloat,3,1>(37./160., 5./24., 1./6.);

   // Set the gamma parameters
   const Eigen::Array<MHDFloat,3,1> ImExRK3::mGamma = Eigen::Array<MHDFloat,3,1>(8./15., 5./12., 3./4.);

   // Set the zeta parameters
   const Eigen::Array<MHDFloat,3,1> ImExRK3::mZeta = Eigen::Array<MHDFloat,3,1>(0., -17./60., -5./12.);

   // Set the zeta parameters
   const Eigen::Array<MHDFloat,3,1> ImExRK3::mCEx = Eigen::Array<MHDFloat,3,1>(0., 0., 1.);

   MHDFloat ImExRK3::lhsT(const int step)
   {
      assert(step < ImExRK3::STEPS);

      return 1./ImExRK3::mBeta(step);
   }

   MHDFloat ImExRK3::lhsL(const int step)
   {
      assert(step < ImExRK3::STEPS);

      return 1.0;
   }

   MHDFloat ImExRK3::rhsT(const int i, const int step)
   {
      assert(step < ImExRK3::STEPS);
      assert(i > -1);
      assert(i < ImExRK3::FIELD_MEMORY+1);

      return 1./ImExRK3::mBeta(step);
   }

   MHDFloat ImExRK3::rhsL(const int i, const int step)
   {
      assert(step < ImExRK3::STEPS);
      assert(i > -1);
      assert(i < ImExRK3::FIELD_MEMORY+1);

      return ImExRK3::mAlpha(step)/ImExRK3::mBeta(step);
   }

   MHDFloat ImExRK3::rhsN(const int i, const int step)
   {
      assert(step < ImExRK3::STEPS);
      assert(i > -1);
      assert(i < ImExRK3::NONLINEAR_MEMORY+1);

      MHDFloat coeff;
      if(i == 0)
      {
         coeff = ImExRK3::mGamma(step)/ImExRK3::mBeta(step);
      } else
      {
         coeff = ImExRK3::mZeta(step)/ImExRK3::mBeta(step);
      }

      return coeff;
   }

   MHDFloat ImExRK3::cEx(const int i)
   {
      return ImExRK3::mCEx(i);
   }

   int ImExRK3::fieldMemory(const int step)
   {
      assert(step < ImExRK3::STEPS);

      return ImExRK3::FIELD_MEMORY;
   }

   int ImExRK3::nonlinearMemory(const int step)
   {
      assert(step < ImExRK3::STEPS);

      int mem;

      // First step does not use t_(n-1) nonlinear term
      if(step == 0)
      {
         mem = 0;
      } else
      {
         mem = ImExRK3::NONLINEAR_MEMORY;
      }

      return mem;
   }

   void ImExRK3::init()
   {
   }

}
}
