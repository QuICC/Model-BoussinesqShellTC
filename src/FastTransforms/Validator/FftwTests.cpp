/**
 * @file FftwTests.cpp
 * @brief Collection of tests to validate FFTW agains known issues
 */

// System includes
//
#include "Debug/DebuggerMacro.h"

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
      fftw_plan plan = fftw_plan_dft_r2c_1d(in.size(),in.data(),reinterpret_cast<fftw_complex*>(out.data()),FftwLibrary::planFlag());

      /* initialise input */
      const MHDFloat dk = static_cast<MHDFloat>(k);
      for(int i = 0;i < in.size();i++)
      {
         MHDFloat g = 2.0*Math::PI*static_cast<MHDFloat>(i)/static_cast<MHDFloat>(in.size());
         in(i) = 2.0*dk*(std::cos(dk*g) + std::sin(dk*g));
      }

      // compute FFT
      fftw_execute(plan);
      fftw_destroy_plan(plan);

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

      DebuggerMacro_showValue("Validate R2C error: ", 1, max_err);
      if(max_err > std::numeric_limits<MHDFloat>::epsilon()*10)
      {
         throw std::logic_error("FFTW R2C validation failed!");
      }
   }

   void FftwTests::validateC2R()
   {
      const int k = 1;
      const int size = 32;
      ArrayZ in(size/2 + 1);
      Array out(size);

      // Create forward plan (R2C)
      fftw_plan plan = fftw_plan_dft_c2r_1d(out.size(),reinterpret_cast<fftw_complex*>(in.data()),out.data(),FftwLibrary::planFlag());

      /* initialise input */
      const MHDFloat dk = static_cast<MHDFloat>(k);
      in.setConstant(0.0);
      in(k) = MHDComplex(1.0, -1.0);

      // compute FFT
      fftw_execute(plan);
      fftw_destroy_plan(plan);

      // Compute error
      MHDFloat max_err = -1.0;
      MHDFloat dsize = static_cast<MHDFloat>(out.size());
      for(int i = 0;i < out.size();i++)
      {
         MHDFloat g = 2.0*Math::PI*static_cast<MHDFloat>(i)/dsize;
         max_err =  std::max(max_err, std::abs(out(i) - 2.0*dk*(std::cos(dk*g) + std::sin(dk*g))));
      }

      DebuggerMacro_showValue("Validate C2R error: ", 1, max_err);
      if(max_err > std::numeric_limits<MHDFloat>::epsilon()*10)
      {
         throw std::logic_error("FFTW C2R validation failed!");
      }
   }

   void FftwTests::validateDftFwd()
   {
      const int k = 1;
      const int size = 32;
      ArrayZ in(size);
      ArrayZ out(size);

      // Create forward plan
      fftw_plan plan = fftw_plan_dft_1d(size, reinterpret_cast<fftw_complex*>(in.data()), reinterpret_cast<fftw_complex*>(out.data()), FFTW_FORWARD, FftwLibrary::planFlag());

      in.setZero();
      out.setZero();

      /* initialise input */
      const MHDFloat dk = static_cast<MHDFloat>(k);
      for(int i = 0;i < in.size();i++)
      {
         MHDFloat g = 2.0*Math::PI*static_cast<MHDFloat>(i)/static_cast<MHDFloat>(in.size());
         in(i) = MHDComplex(3.0, -8.0)*std::exp(Math::cI*dk*g);
         in(i) += MHDComplex(-5.5, 7.1)*std::exp(-Math::cI*dk*g);
      }

      // compute FFT
      fftw_execute(plan);
      fftw_destroy_plan(plan);

      // Compute error
      MHDFloat max_err = -1.0;
      MHDFloat dsize = static_cast<MHDFloat>(size);
      for(int i = 0;i < out.size();i++)
      {
         MHDFloat err = 1.0;
         if(i == k)
         {
            err = std::abs(MHDComplex(3.0, -8.0) - out(i)/dsize);
         } else if(i > 0 && k == size - i)
         {
            err = std::abs(MHDComplex(-5.5, 7.1) - out(i)/dsize);
         } else
         {
            err = std::abs(out(i)/dsize);
         }
         max_err = std::max(max_err, std::abs(err));
      }

      DebuggerMacro_showValue("Validate DFT forward error: ", 1, max_err);
      if(max_err > std::numeric_limits<MHDFloat>::epsilon()*10)
      {
         throw std::logic_error("FFTW DFT forward validation failed!");
      }
   }

   void FftwTests::validateDftBwd()
   {
      const int k = 1;
      const int size = 32;
      ArrayZ in(size);
      ArrayZ out(size);

      // Create forward plan
      fftw_plan plan = fftw_plan_dft_1d(size, reinterpret_cast<fftw_complex*>(in.data()), reinterpret_cast<fftw_complex*>(out.data()), FFTW_BACKWARD, FftwLibrary::planFlag());

      in.setZero();
      out.setZero();

      /* initialise input */
      in(k) = MHDComplex(3.0, -8.0);
      if(k > 0)
      {
         in(size-k) = MHDComplex(-5.5, 7.1);
      }

      // compute FFT
      fftw_execute(plan);
      fftw_destroy_plan(plan);

      // Compute error
      MHDFloat max_err = -1.0;
      const MHDFloat dk = static_cast<MHDFloat>(k);
      for(int i = 0;i < out.size();i++)
      {
         MHDFloat g = 2.0*Math::PI*static_cast<MHDFloat>(i)/static_cast<MHDFloat>(in.size());
         MHDComplex ref = MHDComplex(3.0, -8.0)*std::exp(Math::cI*dk*g);
         ref += MHDComplex(-5.5, 7.1)*std::exp(-Math::cI*dk*g);

         max_err = std::max(max_err, std::abs(std::abs(out(i) - ref)/std::abs(ref)));
      }

      DebuggerMacro_showValue("Validate DFT backward error: ", 1, max_err);
      DebuggerMacro_showValue("Validate DFT backward max OK error: ", 1, std::numeric_limits<MHDFloat>::epsilon()*10);
      if(max_err > std::numeric_limits<MHDFloat>::epsilon()*1e2)
      {
         throw std::logic_error("FFTW DFT backward validation failed!");
      }
   }

   void FftwTests::validateR2R10()
   {
      // Compute error
      MHDFloat max_err = -1.0;

      DebuggerMacro_showValue("Validate R2R REDFT10 error: ", 1, max_err);
      if(max_err > std::numeric_limits<MHDFloat>::epsilon()*10)
      {
         throw std::logic_error("FFTW R2R REDFT10 validation failed!");
      }
   }

   void FftwTests::validateR2R01()
   {
      // Compute error
      MHDFloat max_err = -1.0;

      DebuggerMacro_showValue("Validate R2R REDFT01 error: ", 1, max_err);
      if(max_err > std::numeric_limits<MHDFloat>::epsilon()*10)
      {
         throw std::logic_error("FFTW R2R REDFT01 validation failed!");
      }
   }

   void FftwTests::validateR2R11()
   {
      // Compute error
      MHDFloat max_err = -1.0;

      DebuggerMacro_showValue("Validate R2R REDFT11 error: ", 1, max_err);
      if(max_err > std::numeric_limits<MHDFloat>::epsilon()*10)
      {
         throw std::logic_error("FFTW R2R REDFT11 validation failed!");
      }
   }

}
}
}
