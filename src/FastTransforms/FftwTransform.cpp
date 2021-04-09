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
#include "Exceptions/Exception.hpp"
#include "Base/MathConstants.hpp"
#include "FastTransforms/FftwLibrary.hpp"
#include "FastTransforms/Validator/FftwTests.hpp"

namespace QuICC {

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

      // Run validation tests
      Validator::FftwTests::validateR2C();
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

      // Extract block sizes of mean values
      if(this->mspSetup->type() == FftSetup::COMPLEX)
      {
         int start = 0;
         for(int i = 0; i < this->mspSetup->idBlocks().rows(); ++i)
         {
            if(this->mspSetup->idBlocks()(i,0) == 0)
            {
               this->mMeanBlocks.push_back(std::make_pair(start, this->mspSetup->idBlocks()(i,1)));
            }

            start += this->mspSetup->idBlocks()(i,1);
         }
      }
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

   void FftwTransform::integrate(MatrixZ& rFFTVal, const Matrix& physVal, FftwTransform::IntegratorType::Id integrator)
   {
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
         ArrayZ factor = this->mspSetup->scale()*this->mspSetup->boxScale()*Math::cI*Array::LinSpaced(this->mspSetup->bwdSize(), 0, this->mspSetup->bwdSize()-1);

         // Rescale results
         rFFTVal = factor.asDiagonal()*rFFTVal;

      // Compute first derivative integration
      } else if(integrator == FftwTransform::IntegratorType::INTGDIFF2)
      {
         // Get differentiation factors
         ArrayZ factor = this->mspSetup->scale()*this->mspSetup->boxScale()*Math::cI*Array::LinSpaced(this->mspSetup->bwdSize(), 0, this->mspSetup->bwdSize()-1);
         factor = factor.array().pow(2);

         // Rescale results
         rFFTVal = factor.asDiagonal()*rFFTVal;
      
      // Compute first derivative integration and mean (k = 0 is not zeroed)
      } else if(integrator == FftwTransform::IntegratorType::INTGDIFFM)
      {
         // Get differentiation factors
         ArrayZ factor = this->mspSetup->scale()*this->mspSetup->boxScale()*Math::cI*Array::LinSpaced(this->mspSetup->bwdSize(), 0, this->mspSetup->bwdSize()-1);
         factor(0) = this->mspSetup->scale();

         // Rescale results
         rFFTVal = factor.asDiagonal()*rFFTVal;
      
      // Compute first derivative integration and negative mean (k = 0 is not zeroed)
      } else if(integrator == FftwTransform::IntegratorType::INTGDIFFNEGM)
      {
         // Get differentiation factors
         ArrayZ factor = this->mspSetup->scale()*this->mspSetup->boxScale()*Math::cI*Array::LinSpaced(this->mspSetup->bwdSize(), 0, this->mspSetup->bwdSize()-1);
         factor(0) = -this->mspSetup->scale();

         // Rescale results
         rFFTVal = factor.asDiagonal()*rFFTVal;

      // Compute simple projection
      } else
      {
         // Rescale output from FFT
         rFFTVal *= this->mspSetup->scale();
      }
   }

   void FftwTransform::project(Matrix& rPhysVal, const MatrixZ& fftVal, FftwTransform::ProjectorType::Id projector)
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

      } else if(projector == FftwTransform::ProjectorType::DIFF3)
      {
         // Get differentiation factors
         ArrayZ factor = this->mspSetup->boxScale()*Math::cI*Array::LinSpaced(this->mspSetup->specSize(), 0, this->mspSetup->specSize()-1);
         factor = factor.array().pow(3);

         // Rescale results
         this->mTmpRIn.topRows(this->mspSetup->specSize()) = factor.asDiagonal()*fftVal.topRows(this->mspSetup->specSize());

      // Compute mean only projection
      } else if(projector == FftwTransform::ProjectorType::PROJMEANONLY)
      {
         throw Exception("Mean only call is not possible here!");

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

      // Do transform
      fftw_execute_dft_c2r(this->mBPlan, reinterpret_cast<fftw_complex* >(this->mTmpRIn.data()), rPhysVal.data());
   }

   void FftwTransform::integrate(MatrixZ& rFFTVal, const MatrixZ& physVal, FftwTransform::IntegratorType::Id integrator)
   {
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

      // Get size of positive and negative frequency parts
      int negN = this->mspSetup->specSize()/2;
      int posN = negN + (this->mspSetup->specSize()%2);

      // Compute first derivative integration
      if(integrator == FftwTransform::IntegratorType::INTGDIFF)
      {
         // Get differentiation factors
         ArrayZ factor = this->mspSetup->scale()*this->mspSetup->boxScale()*Math::cI*Array::LinSpaced(posN, 0, posN-1);
         ArrayZ rfactor = this->mspSetup->scale()*this->mspSetup->boxScale()*Math::cI*(Array::LinSpaced(negN, 0, negN-1).array() - static_cast<MHDFloat>(negN));

         // Split positive and negative frequencies and compute derivative
         rFFTVal.topRows(posN) = factor.asDiagonal()*rFFTVal.topRows(posN);
         rFFTVal.bottomRows(negN) = rfactor.asDiagonal()*rFFTVal.bottomRows(negN);

      // Compute second derivative integration
      } else if(integrator == FftwTransform::IntegratorType::INTGDIFF2)
      {
         // Get differentiation factors
         Array factor = -this->mspSetup->scale()*(this->mspSetup->boxScale()*Array::LinSpaced(posN, 0, posN-1)).array().pow(2);
         Array rfactor = -this->mspSetup->scale()*(this->mspSetup->boxScale()*(Array::LinSpaced(negN, 0, negN-1).array() - static_cast<MHDFloat>(negN))).array().pow(2);

         // Split positive and negative frequencies and compute derivative
         rFFTVal.topRows(posN) = factor.asDiagonal()*rFFTVal.topRows(posN);
         rFFTVal.bottomRows(negN) = rfactor.asDiagonal()*rFFTVal.bottomRows(negN);

      // Compute horizontal laplacian integration
      } else if(integrator == FftwTransform::IntegratorType::INTGLAPLH)
      {
         // Get k1^2 factors
         Array factor = (this->mspSetup->boxScale()*Array::LinSpaced(posN, 0, posN-1)).array().pow(2);
         Array rfactor = (this->mspSetup->boxScale()*(Array::LinSpaced(negN, 0, negN-1).array() - static_cast<MHDFloat>(negN))).array().pow(2);

         int start = 0;
         int negRow = rFFTVal.rows() - negN;
         for(int i = 0; i < this->mspSetup->idBlocks().rows(); ++i)
         {
            MHDFloat k2 = this->mspSetup->idBlocks()(i,0)*this->mspSetup->boxScale();
            k2 *= k2;

            // Split positive and negative frequencies to compute rescaling
            Array factor2 = -this->mspSetup->scale()*(k2 + factor.array());
            rFFTVal.block(0, start, posN, this->mspSetup->idBlocks()(i,1)) = factor2.asDiagonal()*rFFTVal.block(0, start, posN, this->mspSetup->idBlocks()(i,1));
            factor2 = -this->mspSetup->scale()*(k2 + rfactor.array());
            rFFTVal.block(negRow, start, negN, this->mspSetup->idBlocks()(i,1)) = factor2.asDiagonal()*rFFTVal.block(negRow, start, negN, this->mspSetup->idBlocks()(i,1));

            // Increment block counter
            start += this->mspSetup->idBlocks()(i,1);
         }

      // Compute first derivative integration and mean (k1 = k2 = 0 mode is not zeroed)
      } else if(integrator == FftwTransform::IntegratorType::INTGDIFFM)
      {
         // Get differentiation factors
         ArrayZ factor = this->mspSetup->scale()*this->mspSetup->boxScale()*Math::cI*Array::LinSpaced(posN, 0, posN-1);
         ArrayZ rfactor = this->mspSetup->scale()*this->mspSetup->boxScale()*Math::cI*(Array::LinSpaced(negN, 0, negN-1).array() - static_cast<MHDFloat>(negN));

         // Extract the mean
         std::vector<ArrayZ> mean;
         for(std::vector<std::pair<int,int> >::const_iterator it = this->mMeanBlocks.begin(); it != this->mMeanBlocks.end(); ++it)
         {
            mean.push_back(rFFTVal.block(0, it->first, 1, it->second).transpose());
         }

         // Split positive and negative frequencies and compute derivative
         rFFTVal.topRows(posN) = factor.asDiagonal()*rFFTVal.topRows(posN);
         rFFTVal.bottomRows(negN) = rfactor.asDiagonal()*rFFTVal.bottomRows(negN);

         // Set the mean
         for(size_t j = 0; j < mean.size(); ++j)
         {
            rFFTVal.block(0, this->mMeanBlocks.at(j).first, 1, this->mMeanBlocks.at(j).second) = this->mspSetup->scale()*mean.at(j).transpose();
         }

      // Compute first derivative integration and mean (k1 = k2 = 0 mode is not zeroed)
      } else if(integrator == FftwTransform::IntegratorType::INTGDIFFNEGM)
      {
         // Get differentiation factors
         ArrayZ factor = this->mspSetup->scale()*this->mspSetup->boxScale()*Math::cI*Array::LinSpaced(posN, 0, posN-1);
         ArrayZ rfactor = this->mspSetup->scale()*this->mspSetup->boxScale()*Math::cI*(Array::LinSpaced(negN, 0, negN-1).array() - static_cast<MHDFloat>(negN));

         // Extract the mean
         std::vector<ArrayZ> mean;
         for(std::vector<std::pair<int,int> >::const_iterator it = this->mMeanBlocks.begin(); it != this->mMeanBlocks.end(); ++it)
         {
            mean.push_back(rFFTVal.block(0, it->first, 1, it->second).transpose());
         }

         // Split positive and negative frequencies and compute derivative
         rFFTVal.topRows(posN) = factor.asDiagonal()*rFFTVal.topRows(posN);
         rFFTVal.bottomRows(negN) = rfactor.asDiagonal()*rFFTVal.bottomRows(negN);

         // Set the mean
         for(size_t j = 0; j < mean.size(); ++j)
         {
            rFFTVal.block(0, this->mMeanBlocks.at(j).first, 1, this->mMeanBlocks.at(j).second) = -this->mspSetup->scale()*mean.at(j).transpose();
         }

      // Compute inverse horizontal gradient integration
      } else if(integrator == FftwTransform::IntegratorType::INTGDIFFFINVLAPLH)
      {
         // Get k1 factors
         ArrayZ factor = this->mspSetup->boxScale()*Math::cI*Array::LinSpaced(posN, 0, posN-1);
         ArrayZ rfactor = this->mspSetup->boxScale()*Math::cI*(Array::LinSpaced(negN, 0, negN-1).array() - static_cast<MHDFloat>(negN));
         // Get k1^2 factors
         Array lfactor = (this->mspSetup->boxScale()*Array::LinSpaced(posN, 0, posN-1)).array().pow(2);
         Array rlfactor = (this->mspSetup->boxScale()*(Array::LinSpaced(negN, 0, negN-1).array() - static_cast<MHDFloat>(negN))).array().pow(2);

         int start = 0;
         int negRow = rFFTVal.rows() - negN;
         for(int i = 0; i < this->mspSetup->idBlocks().rows(); ++i)
         {
            MHDFloat k2 = this->mspSetup->idBlocks()(i,0)*this->mspSetup->boxScale();
            k2 *= k2;

            ArrayZ factor2 = -this->mspSetup->scale()*factor.array()*(k2 + lfactor.array()).inverse().array();
            // Fix k1 == k2 == 0 mode
            if(this->mspSetup->idBlocks()(i,0) == 0)
            {
               factor2(0) = 0.0;
            }
            // Split positive and negative frequencies to compute rescaling
            rFFTVal.block(0, start, posN, this->mspSetup->idBlocks()(i,1)) = factor2.asDiagonal()*rFFTVal.block(0, start, posN, this->mspSetup->idBlocks()(i,1));
            factor2 = -this->mspSetup->scale()*rfactor.array()*(k2 + rlfactor.array()).inverse();
            rFFTVal.block(negRow, start, negN, this->mspSetup->idBlocks()(i,1)) = factor2.asDiagonal()*rFFTVal.block(negRow, start, negN, this->mspSetup->idBlocks()(i,1));

            // Increment block counter
            start += this->mspSetup->idBlocks()(i,1);
         }

      // Compute inverse horizontal laplacian integration
      } else if(integrator == FftwTransform::IntegratorType::INTGINVLAPLH)
      {
         // Get k1^2 factors
         Array factor = (this->mspSetup->boxScale()*Array::LinSpaced(posN, 0, posN-1)).array().pow(2);
         Array rfactor = (this->mspSetup->boxScale()*(Array::LinSpaced(negN, 0, negN-1).array() - static_cast<MHDFloat>(negN))).array().pow(2);

         int start = 0;
         int negRow = rFFTVal.rows() - negN;
         for(int i = 0; i < this->mspSetup->idBlocks().rows(); ++i)
         {
            MHDFloat k2 = this->mspSetup->idBlocks()(i,0)*this->mspSetup->boxScale();
            k2 *= k2;

            // Split positive and negative frequencies to compute rescaling
            Array factor2 = -this->mspSetup->scale()*(k2 + factor.array()).inverse();
            // Fix k1 == k2 == 0 mode
            if(this->mspSetup->idBlocks()(i,0) == 0)
            {
               factor2(0) = 0.0;
            }
            rFFTVal.block(0, start, posN, this->mspSetup->idBlocks()(i,1)) = factor2.asDiagonal()*rFFTVal.block(0, start, posN, this->mspSetup->idBlocks()(i,1));
            factor2 = -this->mspSetup->scale()*(k2 + rfactor.array()).inverse();
            rFFTVal.block(negRow, start, negN, this->mspSetup->idBlocks()(i,1)) = factor2.asDiagonal()*rFFTVal.block(negRow, start, negN, this->mspSetup->idBlocks()(i,1));

            // Increment block counter
            start += this->mspSetup->idBlocks()(i,1);
         }

      // Compute integration and zero k2 = 0, k1 != 0
      } else if(integrator == FftwTransform::IntegratorType::INTGM)
      {
         // Zero the mean
         for(std::vector<std::pair<int,int> >::const_iterator it = this->mMeanBlocks.begin(); it != this->mMeanBlocks.end(); ++it)
         {
            rFFTVal.block(1, it->first, rFFTVal.rows()-1, it->second).setZero();
         }

         // Rescale output from FFT
         rFFTVal *= this->mspSetup->scale();

      // Compute integration and zero k2 = 0, k1 != 0
      } else if(integrator == FftwTransform::IntegratorType::INTGMEANONLY)
      {
         // Z
         int start = 0;
         int rows = rFFTVal.rows();
         for(int i = 0; i < this->mspSetup->idBlocks().rows(); ++i)
         {
            if(this->mspSetup->idBlocks()(i,0) == 0)
            {
               rFFTVal.block(0, start, 1, this->mspSetup->idBlocks()(i,1)) *= this->mspSetup->scale();
               rFFTVal.block(1, start, rows-1, this->mspSetup->idBlocks()(i,1)).setZero();
            } else
            {
               rFFTVal.block(0, start, rows, this->mspSetup->idBlocks()(i,1)).setZero();
            }

            // Increment block counter
            start += this->mspSetup->idBlocks()(i,1);
         }

      // Compute simple projection
      } else
      {
         // Rescale output from FFT
         rFFTVal *= this->mspSetup->scale();
      }
   }

   void FftwTransform::project(MatrixZ& rPhysVal, const MatrixZ& fftVal, FftwTransform::ProjectorType::Id projector)
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

      // Compute second derivative
      } else if(projector == FftwTransform::ProjectorType::DIFF2)
      {
         // Get differentiation factors
         Array factor = -(this->mspSetup->boxScale()*Array::LinSpaced(posN, 0, posN-1)).array().pow(2);
         Array rfactor = -(this->mspSetup->boxScale()*(Array::LinSpaced(negN, 0, negN-1).array() - static_cast<MHDFloat>(negN))).array().pow(2);

         // Split positive and negative frequencies and compute derivative
         this->mTmpZIn.topRows(posN) = factor.asDiagonal()*fftVal.topRows(posN);
         this->mTmpZIn.bottomRows(negN) = rfactor.asDiagonal()*fftVal.bottomRows(negN);

      // Compute second derivative
      } else if(projector == FftwTransform::ProjectorType::DIFF3)
      {
         // Get differentiation factors
         Array factor = -(this->mspSetup->boxScale()*Array::LinSpaced(posN, 0, posN-1)).array().pow(3);
         Array rfactor = -(this->mspSetup->boxScale()*(Array::LinSpaced(negN, 0, negN-1).array() - static_cast<MHDFloat>(negN))).array().pow(3);

         // Split positive and negative frequencies and compute derivative
         this->mTmpZIn.topRows(posN) = factor.asDiagonal()*fftVal.topRows(posN);
         this->mTmpZIn.bottomRows(negN) = rfactor.asDiagonal()*fftVal.bottomRows(negN);

      // Compute mean only projection
      } else if(projector == FftwTransform::ProjectorType::PROJMEANONLY)
      {
         // Initialize to zero
         this->mTmpZIn.setZero();

         // Set the mean
         std::vector<ArrayZ> mean;
         for(std::vector<std::pair<int,int> >::const_iterator it = this->mMeanBlocks.begin(); it != this->mMeanBlocks.end(); ++it)
         {
            this->mTmpZIn.block(0, it->first, 1, it->second) = fftVal.block(0, it->first, 1, it->second);
         }

      // Compute horizontal Laplacian
      } else if(projector == FftwTransform::ProjectorType::LAPLH)
      {
         // Get k1^2 factors
         Array factor = (this->mspSetup->boxScale()*Array::LinSpaced(posN, 0, posN-1)).array().pow(2);
         Array rfactor = (this->mspSetup->boxScale()*(Array::LinSpaced(negN, 0, negN-1).array() - static_cast<MHDFloat>(negN))).array().pow(2);

         int start = 0;
         int negRow = fftVal.rows() - negN;
         for(int i = 0; i < this->mspSetup->idBlocks().rows(); ++i)
         {
            MHDFloat k2 = this->mspSetup->idBlocks()(i,0)*this->mspSetup->boxScale();
            k2 *= k2;

            // Split positive and negative frequencies to compute rescaling
            Array factor2 = -(k2 + factor.array());
            this->mTmpZIn.block(0, start, posN, this->mspSetup->idBlocks()(i,1)) = factor2.asDiagonal()*fftVal.block(0, start, posN, this->mspSetup->idBlocks()(i,1));
            factor2 = -(k2 + rfactor.array());
            this->mTmpZIn.block(negRow, start, negN, this->mspSetup->idBlocks()(i,1)) = factor2.asDiagonal()*fftVal.block(negRow, start, negN, this->mspSetup->idBlocks()(i,1));

            // Increment block counter
            start += this->mspSetup->idBlocks()(i,1);
         }

      // Compute fast derivative of horizontal Laplacian
      } else if(projector == FftwTransform::ProjectorType::DIFFFLAPLH)
      {
         // Get k1 factors
         ArrayZ factor1 = (this->mspSetup->boxScale()*Math::cI*Array::LinSpaced(posN, 0, posN-1));
         ArrayZ rfactor1 = (this->mspSetup->boxScale()*Math::cI*(Array::LinSpaced(negN, 0, negN-1).array() - static_cast<MHDFloat>(negN)));

         // Get k1^3 factors
         ArrayZ factor3 = (this->mspSetup->boxScale()*Math::cI*Array::LinSpaced(posN, 0, posN-1)).array().pow(3);
         ArrayZ rfactor3 = (this->mspSetup->boxScale()*Math::cI*(Array::LinSpaced(negN, 0, negN-1).array() - static_cast<MHDFloat>(negN))).array().pow(3);

         int start = 0;
         int negRow = fftVal.rows() - negN;
         for(int i = 0; i < this->mspSetup->idBlocks().rows(); ++i)
         {
            MHDFloat k2 = this->mspSetup->idBlocks()(i,0)*this->mspSetup->boxScale();
            k2 *= k2;

            ArrayZ factor2 = (-k2*factor1.array() + factor3.array());
            // Split positive and negative frequencies to compute rescaling
            this->mTmpZIn.block(0, start, posN, this->mspSetup->idBlocks()(i,1)) = factor2.asDiagonal()*fftVal.block(0, start, posN, this->mspSetup->idBlocks()(i,1));
            factor2 = (-k2*rfactor1.array() + rfactor3.array());
            this->mTmpZIn.block(negRow, start, negN, this->mspSetup->idBlocks()(i,1)) = factor2.asDiagonal()*fftVal.block(negRow, start, negN, this->mspSetup->idBlocks()(i,1));

            // Increment block counter
            start += this->mspSetup->idBlocks()(i,1);
         }

      // Compute slow derivative of horizontal Laplacian
      } else if(projector == FftwTransform::ProjectorType::DIFFSLAPLH)
      {
         // Get k1^2 factors
         Array factor = (this->mspSetup->boxScale()*Array::LinSpaced(posN, 0, posN-1)).array().pow(2);
         Array rfactor = (this->mspSetup->boxScale()*(Array::LinSpaced(negN, 0, negN-1).array() - static_cast<MHDFloat>(negN))).array().pow(2);

         int start = 0;
         int negRow = fftVal.rows() - negN;
         for(int i = 0; i < this->mspSetup->idBlocks().rows(); ++i)
         {
            MHDComplex k2 = this->mspSetup->idBlocks()(i,0)*this->mspSetup->boxScale()*Math::cI;

            // Split positive and negative frequencies to compute rescaling
            ArrayZ factor2 = k2*(k2*k2 - factor.cast<MHDComplex>().array());
            this->mTmpZIn.block(0, start, posN, this->mspSetup->idBlocks()(i,1)) = factor2.asDiagonal()*fftVal.block(0, start, posN, this->mspSetup->idBlocks()(i,1));
            factor2 = k2*(k2*k2 - rfactor.cast<MHDComplex>().array());
            this->mTmpZIn.block(negRow, start, negN, this->mspSetup->idBlocks()(i,1)) = factor2.asDiagonal()*fftVal.block(negRow, start, negN, this->mspSetup->idBlocks()(i,1));

            // Increment block counter
            start += this->mspSetup->idBlocks()(i,1);
         }

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

      // Do transform
      fftw_execute_dft(this->mBPlan, reinterpret_cast<fftw_complex *>(this->mTmpZIn.data()), reinterpret_cast<fftw_complex *>(rPhysVal.data()));
   }

   void FftwTransform::forceConjugate(MatrixZ& rFFTVal)
   {
      // Get size of positive and negative frequency parts
      int negN = this->mspSetup->specSize()/2;
      int posN = negN + (this->mspSetup->specSize()%2);
      int endN = rFFTVal.rows();

      // Loop over special blocks
      for(std::vector<std::pair<int,int> >::const_iterator it = this->mMeanBlocks.begin(); it != this->mMeanBlocks.end(); ++it)
      {
         // Copy complex conjugate into negative frequency part
         for(int i = 1; i < posN; i++)
         {
            rFFTVal.block(endN - i, it->first, 1, it->second) = rFFTVal.block(i, it->first, 1, it->second).conjugate();
         }
      }
   }

#ifdef QUICC_STORAGEPROFILE
   MHDFloat FftwTransform::requiredStorage() const
   {
      MHDFloat mem = 0.0;

      // Storage required for the fftw plans 
      mem += 8.0*2.0;

      return mem;
   }
#endif // QUICC_STORAGEPROFILE

}
}
