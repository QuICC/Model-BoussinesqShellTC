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
      if(max_err > std::numeric_limits<MHDFloat>::epsilon()*1e2)
      {
         throw std::logic_error("FFTW DFT backward validation failed!");
      }
   }

   void FftwTests::validateR2R10()
   {
      const int size = 32;
      Array in(size);
      Array out(size);

      // Create forward plan
      fftw_plan plan = fftw_plan_r2r_1d(size, in.data(), out.data(), FFTW_REDFT10, FftwLibrary::planFlag());

      in.setZero();
      out.setZero();

      Array ref(5);
      ref(0) = 0.7;
      ref(1) = -0.35;
      ref(2) = 0.1;
      ref(3) = 0.44;
      ref(4) = -0.27;

      /* initialise input */
      for(int i = 0;i < in.size();i++)
      {
         MHDFloat x = std::cos((Math::PI)*(static_cast<MHDFloat>(i)+0.5)/static_cast<MHDFloat>(size));
         in(i)  = ref(0);
         in(i)  += 2.0*ref(1)*x;
         in(i)  += 2.0*ref(2)*(-1.0 + 2.0*x*x);
         in(i)  += 2.0*ref(3)*(-3.0*x + 4.0*std::pow(x,3));
         in(i)  += 2.0*ref(4)*(1.0 - 8.0*x*x + 8.0*std::pow(x,4));
      }

      // compute FFT
      fftw_execute(plan);
      fftw_destroy_plan(plan);

      // Compute error
      MHDFloat max_err = -1.0;
      MHDFloat dsize = static_cast<MHDFloat>(2*size);
      for(int i = 0;i < out.size();i++)
      {
         if(i < ref.size())
         {
            max_err = std::max(max_err, std::abs(out(i)/dsize - ref(i)));
         } else
         {
            max_err = std::max(max_err, std::abs(out(i)/dsize));
         }
      }

      DebuggerMacro_showValue("Validate R2R REDFT10 error: ", 1, max_err);
      if(max_err > std::numeric_limits<MHDFloat>::epsilon()*10)
      {
         throw std::logic_error("FFTW R2R REDFT10 validation failed!");
      }
   }

   void FftwTests::validateR2R01()
   {
      const int size = 32;
      Array in(size);
      Array out(size);

      // Create forward plan
      fftw_plan plan = fftw_plan_r2r_1d(size, in.data(), out.data(), FFTW_REDFT01, FftwLibrary::planFlag());

      in.setZero();
      out.setZero();

      Array ref(5);
      ref(0) = 0.7;
      ref(1) = -0.35;
      ref(2) = 0.1;
      ref(3) = 0.44;
      ref(4) = -0.27;

      /* initialise input */
      in(0) = ref(0);
      in(1) = ref(1);
      in(2) = ref(2);
      in(3) = ref(3);
      in(4) = ref(4);

      // compute FFT
      fftw_execute(plan);
      fftw_destroy_plan(plan);

      // Compute error
      MHDFloat max_err = -1.0;
      for(int i = 0;i < out.size();i++)
      {
         MHDFloat x = std::cos((Math::PI)*(static_cast<MHDFloat>(i)+0.5)/static_cast<MHDFloat>(size));
         MHDFloat r  = ref(0);
         r  += 2.0*ref(1)*x;
         r  += 2.0*ref(2)*(-1.0 + 2.0*x*x);
         r  += 2.0*ref(3)*(-3.0*x + 4.0*std::pow(x,3));
         r  += 2.0*ref(4)*(1.0 - 8.0*x*x + 8.0*std::pow(x,4));

         max_err = std::max(max_err, std::abs(out(i) - r));
      }

      DebuggerMacro_showValue("Validate R2R REDFT01 error: ", 1, max_err);
      if(max_err > std::numeric_limits<MHDFloat>::epsilon()*10)
      {
         throw std::logic_error("FFTW R2R REDFT01 validation failed!");
      }
   }

   void FftwTests::validateR2R11()
   {
      const int size = 32;
      Array in(size);
      Array out(size);

      // Create forward plan
      fftw_plan plan = fftw_plan_r2r_1d(size, in.data(), out.data(), FFTW_REDFT11, FftwLibrary::planFlag());

      in.setZero();
      out.setZero();

      Array ref(5);
      ref(0) = -0.35;
      ref(1) = 0.1;
      ref(2) = 0.44;
      ref(3) = -0.27;
      ref(4) = -0.13;

      /* initialise input */
      for(int i = 0;i < in.size();i++)
      {
         MHDFloat x = std::cos((Math::PI)*(static_cast<MHDFloat>(i)+0.5)/static_cast<MHDFloat>(size));
         MHDFloat r = std::sqrt((x+1.0)/2.0);
         in(i)  += 2.0*ref(0)*(r);
         in(i)  += 2.0*ref(1)*(r*(-3.0 + 4.0*std::pow(r,2)));
         in(i)  += 2.0*ref(2)*(r*(5.0 - 20.0*std::pow(r,2) + 16.0*std::pow(r,4)));
         in(i)  += 2.0*ref(3)*(r*(-7.0 + 8.0*std::pow(r,2)*(7.0 - 14.0*std::pow(r,2) + 8.0*std::pow(r,4))));
         in(i)  += 2.0*ref(4)*(r*(-3.0 + 4.0*std::pow(r,2))*(-3.0 + 4.0*std::pow(r,2)*std::pow(3.0 - 4.0*std::pow(r,2),2)));
      }

      // compute FFT
      fftw_execute(plan);
      fftw_destroy_plan(plan);

      // Compute error
      MHDFloat max_err = -1.0;
      MHDFloat dsize = static_cast<MHDFloat>(2*size);
      for(int i = 0;i < out.size();i++)
      {
         if(i < ref.size())
         {
            max_err = std::max(max_err, std::abs(out(i)/dsize - ref(i)));
         } else
         {
            max_err = std::max(max_err, std::abs(out(i)/dsize));
         }
      }

      DebuggerMacro_showValue("Validate R2R REDFT11 error: ", 1, max_err);
      if(max_err > std::numeric_limits<MHDFloat>::epsilon()*10)
      {
         throw std::logic_error("FFTW R2R REDFT11 validation failed!");
      }
   }

}
}
}
