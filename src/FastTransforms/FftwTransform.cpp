/** 
 * @file FftwTransform.cpp
 * @brief Source of the implementation of the FFTW transform
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "FastTransforms/FftwTransform.hpp"

// Project includes
//
#include "Base/MathConstants.hpp"
#include "FastTransforms/FftwLibrary.hpp"

namespace GeoMHDiSCC {

namespace Transform {

   Array FftwTransform::generateGrid(const int size)
   {
      // Initialise grid storage
      Array grid(size);

      // Create equispaced FFT grid
      for(int k = 0; k < size; k++)
      {
         grid(k) = 2.0*Math::PI*static_cast<MHDFloat>(k)/static_cast<MHDFloat>(size);
      }

      return grid;
   }

   FftwTransform::FftwTransform()
      : mFPlan(NULL), mBPlan(NULL)
   {
      // Initialize FFTW
      FftwLibrary::initFft();
   }

   FftwTransform::~FftwTransform()
   {
      // Cleanup memory used by FFTW
      FftwLibrary::cleanupFft();
   }

   void FftwTransform::init(FftwTransform::SharedSetupType spSetup)
   {
      // Store the shared pointer to setup object
      this->mspSetup = spSetup;

      // Set the scaling factor
      this->mspSetup->setScale(1.0/static_cast<MHDFloat>(this->mspSetup->fwdSize()));

      // Initialise FFTW plans
      this->initFft();

      // Register the FFTW object
      FftwLibrary::registerFft();
   }

   void FftwTransform::requiredOptions(std::set<NonDimensional::Id>& list, const Dimensions::Transform::Id dimId) const
   {
      //
      // No possible options
      //
   }

   void FftwTransform::setOptions(const std::map<NonDimensional::Id, MHDFloat>& options, const Dimensions::Transform::Id dimId)
   {
      //
      // No possible options
      //
   }

   Array FftwTransform::meshGrid() const
   {
      return FftwTransform::generateGrid(this->mspSetup->fwdSize());
   }

   void FftwTransform::initFft()
   {
      int fwdSize = this->mspSetup->fwdSize();
      int bwdSize = this->mspSetup->bwdSize();
      int howmany = this->mspSetup->howmany();

      // Create the two plans
      const int  *fftSize = &fwdSize;

      if(this->mspSetup->type() == FftSetup::MIXED)
      {
         // create temporary storage for plan computation
         Matrix    tmpReal(fwdSize, howmany);
         MatrixZ   tmpCplx(bwdSize, howmany);

         // Create the real to complex plan
         this->mFPlan = fftw_plan_many_dft_r2c(1, fftSize, howmany, tmpReal.data(), NULL, 1, fwdSize, reinterpret_cast<fftw_complex* >(tmpCplx.data()), NULL, 1, bwdSize, FftwLibrary::planFlag());

         // Create the complex to real plan
         this->mBPlan = fftw_plan_many_dft_c2r(1, fftSize, howmany, reinterpret_cast<fftw_complex* >(tmpCplx.data()), NULL, 1, bwdSize, tmpReal.data(), NULL, 1, fwdSize, FftwLibrary::planFlag());
      } else
      {
         MatrixZ   tmpCplxA(fwdSize, howmany);
         MatrixZ   tmpCplxB(bwdSize, howmany);

         // Create the forward complex to complex plan
         this->mFPlan = fftw_plan_many_dft(1, fftSize, howmany, reinterpret_cast<fftw_complex* >(tmpCplxA.data()), NULL, 1, fwdSize, reinterpret_cast<fftw_complex* >(tmpCplxB.data()), NULL, 1, bwdSize, FFTW_FORWARD, FftwLibrary::planFlag());

         // Create the backward complex to complex plan
         this->mBPlan = fftw_plan_many_dft(1, fftSize, howmany, reinterpret_cast<fftw_complex* >(tmpCplxB.data()), NULL, 1, bwdSize, reinterpret_cast<fftw_complex* >(tmpCplxA.data()), NULL, 1, fwdSize, FFTW_BACKWARD, FftwLibrary::planFlag());
      }

      // Initialise temporary storage
      this->mTmpRIn.setZero(bwdSize, howmany);
      this->mTmpROut.setZero(fwdSize, howmany);
      this->mTmpZIn.setZero(bwdSize, howmany);
      this->mTmpROut.setZero(fwdSize, howmany);
   }

   void FftwTransform::cleanupFft()
   {
      // Detroy forward plan
      if(this->mFPlan)
      {
         fftw_destroy_plan(this->mFPlan);
      }

      // Detroy backward plan
      if(this->mBPlan)
      {
         fftw_destroy_plan(this->mBPlan);
      }

      // Unregister the FFTW object
      FftwLibrary::unregisterFft();

      // cleanup FFTW library
      FftwLibrary::cleanupFft();
   }

   void FftwTransform::integrate(MatrixZ& rFFTVal, const Matrix& physVal, FftwTransform::IntegratorType::Id integrator, Arithmetics::Id arithId)
   {
      assert(arithId == Arithmetics::SET);

      // Assert that a mixed transform was setup
      assert(this->mspSetup->type() == FftSetup::MIXED);

      // assert right sizes for input matrix
      assert(physVal.rows() == this->mspSetup->fwdSize());
      assert(physVal.cols() == this->mspSetup->howmany());

      // assert right sizes for output matrix
      assert(rFFTVal.rows() == this->mspSetup->bwdSize());
      assert(rFFTVal.cols() == this->mspSetup->howmany());

      // Do transform
      fftw_execute_dft_r2c(this->mFPlan, const_cast<MHDFloat *>(physVal.data()), reinterpret_cast<fftw_complex* >(rFFTVal.data()));

      // Compute first derivative integration
      if(integrator == FftwTransform::IntegratorType::INTGDIFF)
      {
         // Get differentiation factors
         ArrayZ factor = -this->mspSetup->scale()*this->mspSetup->boxScale()*Math::cI*Array::LinSpaced(this->mspSetup->bwdSize(), 0, this->mspSetup->bwdSize()-1);

         // Rescale results
         rFFTVal = factor.asDiagonal()*rFFTVal;

      // Compute simple projection
      } else
      {
         // Rescale output from FFT
         rFFTVal *= this->mspSetup->scale();
      }
   }

   void FftwTransform::project(Matrix& rPhysVal, const MatrixZ& fftVal, FftwTransform::ProjectorType::Id projector, Arithmetics::Id arithId)
   {
      // Assert that a mixed transform was setup
      assert(this->mspSetup->type() == FftSetup::MIXED);

      // assert on the padding size
      assert(this->mspSetup->padSize() >= 0);
      assert(this->mspSetup->bwdSize() - this->mspSetup->padSize() >= 0);

      // assert right sizes for input  matrix
      assert(fftVal.rows() == this->mspSetup->bwdSize());
      assert(fftVal.cols() == this->mspSetup->howmany());

      // assert right sizes for output matrix
      assert(rPhysVal.rows() == this->mspSetup->fwdSize());
      assert(rPhysVal.cols() == this->mspSetup->howmany());

      // Compute first derivative
      if(projector == FftwTransform::ProjectorType::DIFF)
      {
         // Get differentiation factors
         ArrayZ factor = this->mspSetup->boxScale()*Math::cI*Array::LinSpaced(this->mspSetup->specSize(), 0, this->mspSetup->specSize()-1);

         // Rescale results
         this->mTmpRIn.topRows(this->mspSetup->specSize()) = factor.asDiagonal()*fftVal.topRows(this->mspSetup->specSize());

      } else if(projector == FftwTransform::ProjectorType::DIFF2)
      {
         // Get differentiation factors
         ArrayZ factor = this->mspSetup->boxScale()*Math::cI*Array::LinSpaced(this->mspSetup->specSize(), 0, this->mspSetup->specSize()-1);
         factor = factor.array().pow(2);

         // Rescale results
         this->mTmpRIn.topRows(this->mspSetup->specSize()) = factor.asDiagonal()*fftVal.topRows(this->mspSetup->specSize());

      // Compute simple projection
      } else
      {
         // Rescale results
         this->mTmpRIn.topRows(this->mspSetup->specSize()) = fftVal.topRows(this->mspSetup->specSize());
      }

      // Set the m=0 values to zero
      this->mTmpRIn.row(0).imag().setConstant(0);

      // Set the padded values to zero
      this->mTmpRIn.bottomRows(this->mspSetup->padSize()).setZero();

      // Do transform with corresponding arithmetics
      if(arithId == Arithmetics::SET)
      {
         fftw_execute_dft_c2r(this->mBPlan, reinterpret_cast<fftw_complex* >(this->mTmpRIn.data()), rPhysVal.data());
      } else if(arithId == Arithmetics::SETNEG)
      {
         this->mTmpRIn = -this->mTmpRIn;
         fftw_execute_dft_c2r(this->mBPlan, reinterpret_cast<fftw_complex* >(this->mTmpRIn.data()), rPhysVal.data());
      } else if(arithId == Arithmetics::ADD)
      {
         fftw_execute_dft_c2r(this->mBPlan, reinterpret_cast<fftw_complex* >(this->mTmpRIn.data()), this->mTmpROut.data());
         rPhysVal += this->mTmpROut;
      } else if(arithId == Arithmetics::SUB)
      {
         fftw_execute_dft_c2r(this->mBPlan, reinterpret_cast<fftw_complex* >(this->mTmpRIn.data()), this->mTmpROut.data());
         rPhysVal -= this->mTmpROut;
      }
   }

   void FftwTransform::integrate(MatrixZ& rFFTVal, const MatrixZ& physVal, FftwTransform::IntegratorType::Id integrator, Arithmetics::Id arithId)
   {
      assert(arithId == Arithmetics::SET);

      // Assert that a non mixed transform was setup
      assert(this->mspSetup->type() == FftSetup::COMPLEX);

      // assert right sizes for input matrix
      assert(physVal.rows() == this->mspSetup->fwdSize());
      assert(physVal.cols() == this->mspSetup->howmany());

      // assert right sizes for output matrix
      assert(rFFTVal.rows() == this->mspSetup->bwdSize());
      assert(rFFTVal.cols() == this->mspSetup->howmany());

      // Do transform
      fftw_execute_dft(this->mFPlan, reinterpret_cast<fftw_complex* >(const_cast<MHDComplex *>(physVal.data())), reinterpret_cast<fftw_complex* >(rFFTVal.data()));

      // Compute first derivative integration
      if(integrator == FftwTransform::IntegratorType::INTGDIFF)
      {
         // Get differentiation factors
         ArrayZ factor = this->mspSetup->scale()*this->mspSetup->boxScale()*Math::cI*Array::LinSpaced(this->mspSetup->bwdSize(), 0, this->mspSetup->bwdSize()-1);

         // Rescale results
         rFFTVal = factor.asDiagonal()*rFFTVal;

      // Compute simple projection
      } else
      {
         // Rescale output from FFT
         rFFTVal *= this->mspSetup->scale();
      }
   }

   void FftwTransform::project(MatrixZ& rPhysVal, const MatrixZ& fftVal, FftwTransform::ProjectorType::Id projector, Arithmetics::Id arithId)
   {
      // Assert that a non mixed transform was setup
      assert(this->mspSetup->type() == FftSetup::COMPLEX);

      // assert on the padding size
      assert(this->mspSetup->padSize() >= 0);
      assert(this->mspSetup->bwdSize() - this->mspSetup->padSize() >= 0);

      // assert right sizes for input  matrix
      assert(fftVal.rows() == this->mspSetup->bwdSize());
      assert(fftVal.cols() == this->mspSetup->howmany());

      // assert right sizes for output matrix
      assert(rPhysVal.rows() == this->mspSetup->fwdSize());
      assert(rPhysVal.cols() == this->mspSetup->howmany());

      // Get size of positive and negative frequency parts
      int negN = this->mspSetup->specSize()/2;
      int posN = negN + (this->mspSetup->specSize()%2);

      // Compute first derivative
      if(projector == FftwTransform::ProjectorType::DIFF)
      {
         // Get differentiation factors
         ArrayZ factor = this->mspSetup->boxScale()*Math::cI*Array::LinSpaced(posN, 0, posN-1);
         ArrayZ rfactor = this->mspSetup->boxScale()*Math::cI*(Array::LinSpaced(negN, 0, negN-1).array() - static_cast<MHDFloat>(negN));

         // Split positive and negative frequencies and compute derivative
         this->mTmpZIn.topRows(posN) = factor.asDiagonal()*fftVal.topRows(posN);
         this->mTmpZIn.bottomRows(negN) = rfactor.asDiagonal()*fftVal.bottomRows(negN);

      // Compute first derivative
      } else if(projector == FftwTransform::ProjectorType::DIFF2)
      {
         // Get differentiation factors
         Array factor = -(this->mspSetup->boxScale()*Array::LinSpaced(posN, 0, posN-1)).array().pow(2);
         Array rfactor = -(this->mspSetup->boxScale()*(Array::LinSpaced(negN, 0, negN-1).array() - static_cast<MHDFloat>(negN))).array().pow(2);

         // Split positive and negative frequencies and compute derivative
         this->mTmpZIn.topRows(posN) = factor.asDiagonal()*fftVal.topRows(posN);
         this->mTmpZIn.bottomRows(negN) = rfactor.asDiagonal()*fftVal.bottomRows(negN);

      // Compute simple projection
      } else
      {
         // Split positive and negative frequencies
         this->mTmpZIn.topRows(posN) = fftVal.topRows(posN);
         this->mTmpZIn.bottomRows(negN) = fftVal.bottomRows(negN);
      }

      // Force complex conjugate symmetry for zero modes
      this->forceConjugate(this->mTmpZIn);

      // Set the padded values to zero
      this->mTmpZIn.block(posN, 0, this->mspSetup->padSize(), this->mTmpZIn.cols()).setZero();

      // Do transform with corresponding arithmetics
      if(arithId == Arithmetics::SET)
      {
         fftw_execute_dft(this->mBPlan, reinterpret_cast<fftw_complex *>(this->mTmpZIn.data()), reinterpret_cast<fftw_complex *>(rPhysVal.data()));

      } else if(arithId == Arithmetics::SETNEG)
      {
         this->mTmpZIn = -this->mTmpZIn;
         fftw_execute_dft(this->mBPlan, reinterpret_cast<fftw_complex *>(this->mTmpZIn.data()), reinterpret_cast<fftw_complex *>(rPhysVal.data()));

      } else if(arithId == Arithmetics::ADD)
      {
         fftw_execute_dft(this->mBPlan, reinterpret_cast<fftw_complex *>(this->mTmpZIn.data()), reinterpret_cast<fftw_complex *>(this->mTmpZOut.data()));
         rPhysVal += this->mTmpZOut;

      } else if(arithId == Arithmetics::SUB)
      {
         fftw_execute_dft(this->mBPlan, reinterpret_cast<fftw_complex *>(this->mTmpZIn.data()), reinterpret_cast<fftw_complex *>(this->mTmpZOut.data()));
         rPhysVal -= this->mTmpZOut;
      }
   }

   void FftwTransform::forceConjugate(MatrixZ& rFFTVal)
   {
      // Get size of positive and negative frequency parts
      int negN = this->mspSetup->specSize()/2;
      int posN = negN + (this->mspSetup->specSize()%2);
      int endN = rFFTVal.rows();

      // Get number of special blocks
      int rows = this->mspSetup->specialBlocks().rows();

      // Loop over special blocks
      for(int j = 0; j < rows; j++)
      {
         // Copy complex conjugate into negative frequency part
         for(int i = 1; i < posN; i++)
         {
            rFFTVal.block(endN - i, this->mspSetup->specialBlocks()(j,0), 1, this->mspSetup->specialBlocks()(j,1)) = rFFTVal.block(i, this->mspSetup->specialBlocks()(j,0), 1, this->mspSetup->specialBlocks()(j,1)).conjugate();
         }
      }
   }

#ifdef GEOMHDISCC_STORAGEPROFILE
   MHDFloat FftwTransform::requiredStorage() const
   {
      MHDFloat mem = 0.0;

      // Storage required for the fftw plans 
      mem += 8.0*2.0;

      return mem;
   }
#endif // GEOMHDISCC_STORAGEPROFILE

}
}
