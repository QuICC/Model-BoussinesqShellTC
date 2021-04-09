/**
 * @file FftwTests.cpp
 * @brief Collection of tests to validate FFTW agains known issues
 */

// System includes
//

// External includes
//
#include <fftw3.h>

// Class include
//
#include "FastTransforms/Validator/FftwTests.hpp"

// Project includes
//
#include "Exceptions/Exception.hpp"
#include "Base/MathConstants.hpp"
#include "FastTransforms/FftwLibrary.hpp"

#include <iostream>
namespace QuICC {

namespace Transform {

namespace Validator {

void FftwTests::validateR2C()
{
   const int k = 1;
   const int size = 32;
   Array in(size);
   ArrayZ out(in.size()/2 + 1);

   // Create forward plan (R2C)
   fftw_plan fwd_plan = fftw_plan_dft_r2c_1d(in.size(),in.data(),reinterpret_cast<fftw_complex*>(out.data()),FftwLibrary::planFlag());

   /* initialise input */
   const MHDFloat dk = static_cast<MHDFloat>(k);
   for(int i = 0;i < in.size();i++) 
   {
      MHDFloat g = 2.0*Math::PI*static_cast<MHDFloat>(i)/static_cast<MHDFloat>(in.size());
      in(i) = 2.0*dk*(std::cos(dk*g) + std::sin(dk*g));
   }

   // compute FFT 
   fftw_execute(fwd_plan); 
   fftw_destroy_plan(fwd_plan);

   // Compute error
   MHDFloat max_err = -1.0;
   MHDFloat dsize = static_cast<MHDFloat>(size);
   for(int i = 0;i < out.size();i++) 
   {
      MHDFloat err = 1.0;
      if(i == k)
      {
         err = std::abs(MHDComplex(1.0, -1.0) - out(i)/dsize);
      } else
      {
         err = std::abs(out(i)/dsize);
      }
      max_err = std::max(max_err, std::abs(err));
   }

   if(max_err > std::numeric_limits<MHDFloat>::epsilon()*10)
   {
      throw std::logic_error("FFTW R2C validation failed!");
   }
}

}
}
}
