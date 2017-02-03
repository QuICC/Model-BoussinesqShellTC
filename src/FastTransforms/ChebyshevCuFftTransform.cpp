/** 
 * @file ChebyshevCuFftTransform.cpp
 * @brief Source of the implementation of the Chebyshev cuFFT transform
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "FastTransforms/ChebyshevCuFftTransform.hpp"

// Project includes
//
#include "Base/MathConstants.hpp"
#include "FastTransforms/CuFftLibrary.hpp"
#include "Python/PythonWrapper.hpp"

namespace QuICC {

namespace Transform {

   Array ChebyshevCuFftTransform::generateGrid(const int size)
   {
      // Initialise grid storage
      Array grid(size);

      // Create Chebyshev grid
      for(int k = 0; k < size; k++)
      {
         grid(k) = std::cos((Math::PI)*(static_cast<MHDFloat>(k)+0.5)/static_cast<MHDFloat>(size));
      }

      return grid;
   }

   ChebyshevCuFftTransform::ChebyshevCuFftTransform()
      : mFPlan(), mBPlan(), mCScale(0.0)
   {
   }

   ChebyshevCuFftTransform::~ChebyshevCuFftTransform()
   {
      // Cleanup memory used by cuFFT
      CuFftLibrary::cleanupFft();
   }

   void ChebyshevCuFftTransform::init(ChebyshevCuFftTransform::SharedSetupType spSetup)
   {
      // Store the shared pointer to setup object
      this->mspSetup = spSetup;

      // Set the scaling factor
      this->mspSetup->setScale(1.0/static_cast<MHDFloat>(2*this->mspSetup->fwdSize()));

      // Register the cuFFT object
      CuFftLibrary::registerFft();

      // Initialise cuFFT interface
      this->initFft();

      // Initialise Chebyshev operator(s)
      this->initOperators();
   }

   void ChebyshevCuFftTransform::requiredOptions(std::set<NonDimensional::Id>& list, const Dimensions::Transform::Id dimId) const
   {
      if(dimId == Dimensions::Transform::TRA1D)
      {
         list.insert(NonDimensional::SCALE1D);

      } else if(dimId == Dimensions::Transform::TRA2D)
      {
         list.insert(NonDimensional::SCALE2D);

      } else if(dimId == Dimensions::Transform::TRA3D)
      {
         list.insert(NonDimensional::SCALE3D);
      }
   }

   void ChebyshevCuFftTransform::setOptions(const std::map<NonDimensional::Id, MHDFloat>& options, const Dimensions::Transform::Id dimId)
   {
      if(dimId == Dimensions::Transform::TRA1D)
      {
         this->mCScale = options.find(NonDimensional::SCALE1D)->second;

      } else if(dimId == Dimensions::Transform::TRA2D)
      {
         this->mCScale = options.find(NonDimensional::SCALE2D)->second;

      } else if(dimId == Dimensions::Transform::TRA3D)
      {
         this->mCScale = options.find(NonDimensional::SCALE3D)->second;
      }
   }

   Array ChebyshevCuFftTransform::meshGrid() const
   {
      return ChebyshevCuFftTransform::generateGrid(this->mspSetup->fwdSize());
   }

   void ChebyshevCuFftTransform::initFft()
   {  
      /// \mhdBug implement strided stranforms for complex <-> complex case if possible

      int fwdSize = 2*this->mspSetup->fwdSize();
      int bwdSize = this->mspSetup->bwdSize()+1;
      int howmany = this->mspSetup->howmany();

      // Create the two plans
      int  *fftSize = &fwdSize;

      if(this->mspSetup->type() == FftSetup::COMPONENT)
      {
      } else
      {
      }

      // Initialise temporary storage
      this->mTmpR.setZero(fwdSize, howmany);
      this->mTmpZ.setZero(bwdSize, howmany);

      // Create the physical to spectral plan
      checkCudaErrors(cufftPlanMany(&this->mFPlan, 1, fftSize, NULL, 1, fwdSize, NULL, 1, bwdSize, CUFFT_D2Z, howmany));

      // Create the spectral to physical plan
      checkCudaErrors(cufftPlanMany(&this->mBPlan, 1, fftSize, NULL, 1, bwdSize, NULL, 1, fwdSize, CUFFT_Z2D, howmany));

      // Allocate common device memory
      checkCudaErrors(cudaMalloc((void **)&this->mDevR, sizeof(cufftDoubleReal)*fwdSize*howmany));
      checkCudaErrors(cudaMalloc((void **)&this->mDevZ, sizeof(cufftDoubleComplex)*bwdSize*howmany));

      // Initialise phases
      this->mPhase = (-Math::cI*2.0*Math::PI*((this->mspSetup->fwdSize()-0.5)/(2.0*this->mspSetup->fwdSize()))*Array::LinSpaced(this->mTmpZ.rows(), 0, this->mTmpZ.rows()-1)).array().exp();
      this->mPhase_1 = 1.0/this->mPhase.array();
   }

   void ChebyshevCuFftTransform::initOperators()
   {
      this->mDiff.resize(this->mspSetup->specSize(),this->mspSetup->specSize());

      // Initialise python wrapper
      PythonWrapper::init();
      PythonWrapper::import("quicc.geometry.cartesian.cartesian_1d");

      // Prepare arguments to d1(...) call
      PyObject *pArgs, *pValue;
      pArgs = PyTuple_New(3);
      // ... get operator size
      pValue = PyLong_FromLong(this->mspSetup->specSize());
      PyTuple_SetItem(pArgs, 0, pValue);
      // ... create boundray condition (none)
      pValue = PyDict_New();
      PyDict_SetItem(pValue, PyLong_FromLong(0), PyLong_FromLong(0));
      PyTuple_SetItem(pArgs, 1, pValue);
      // ... set coefficient to 1.0
      pValue = PyFloat_FromDouble(this->mCScale);
      PyTuple_SetItem(pArgs, 2, pValue);

      // Call d1
      PythonWrapper::setFunction("d1");
      pValue = PythonWrapper::callFunction(pArgs);

      // Fill matrix and clenup
      PythonWrapper::fillMatrix(this->mDiff, pValue);
      Py_DECREF(pValue);
      PythonWrapper::finalize();
   }

   void ChebyshevCuFftTransform::cleanupFft()
   {
      // Detroy forward plan
      if(this->mFPlan)
      {
         checkCudaErrors(cufftDestroy(this->mFPlan));
      }

      // Detroy backward plan
      if(this->mBPlan)
      {
         checkCudaErrors(cufftDestroy(this->mBPlan));
      }

      // Free device memory
      cudaFree(this->mDevR);
      cudaFree(this->mDevZ);

      // Unregister the cuFFT object
      CuFftLibrary::unregisterFft();

      // cleanup cuFFT library
      CuFftLibrary::cleanupFft();
   }

   void ChebyshevCuFftTransform::integrate(Matrix& rChebVal, const Matrix& physVal, ChebyshevCuFftTransform::IntegratorType::Id integrator)
   {
      // Assert that a mixed transform was not setup
      assert(this->mspSetup->type() == FftSetup::REAL);

      // assert right sizes for input matrix
      assert(physVal.rows() == this->mspSetup->fwdSize());
      assert(physVal.cols() == this->mspSetup->howmany());

      // assert right sizes for output matrix
      assert(rChebVal.rows() == this->mspSetup->bwdSize());
      assert(rChebVal.cols() == this->mspSetup->howmany());

      // Copy data to perform 2N R2C
      this->mTmpR.topRows(this->mspSetup->fwdSize()) = physVal;
      this->mTmpR.bottomRows(this->mspSetup->fwdSize()) = physVal.colwise().reverse();

      // Do transform
      checkCudaErrors(cudaMemcpy(this->mDevR, this->mTmpR.data(), sizeof(cufftDoubleReal)*this->mTmpR.size(), cudaMemcpyHostToDevice));
      checkCudaErrors(cufftExecD2Z(this->mFPlan, this->mDevR, this->mDevZ));
      checkCudaErrors(cudaMemcpy(this->mTmpZ.data(), this->mDevZ, sizeof(cufftDoubleComplex)*this->mTmpZ.size(), cudaMemcpyDeviceToHost));

      // Extract real part and rescale to remove FFT scaling
      this->mTmpZ = this->mPhase_1.asDiagonal()*this->mTmpZ;
      rChebVal = this->mspSetup->scale()*this->mTmpZ.topRows(this->mspSetup->bwdSize()).real();
   }

   void ChebyshevCuFftTransform::project(Matrix& rPhysVal, const Matrix& chebVal, ChebyshevCuFftTransform::ProjectorType::Id projector)
   {
      // Assert that a mixed transform was not setup
      assert(this->mspSetup->type() == FftSetup::REAL);

      // assert on the padding size
      assert(this->mspSetup->padSize() >= 0);
      assert(this->mspSetup->bwdSize() - this->mspSetup->padSize() >= 0);

      // assert right sizes for input  matrix
      assert(chebVal.rows() == this->mspSetup->bwdSize());
      assert(chebVal.cols() == this->mspSetup->howmany());

      // assert right sizes for output matrix
      assert(rPhysVal.rows() == this->mspSetup->fwdSize());
      assert(rPhysVal.cols() == this->mspSetup->howmany());

      // Compute first derivative
      if(projector == ChebyshevCuFftTransform::ProjectorType::DIFF)
      {
         this->mTmpZ.topRows(this->mspSetup->specSize()).real() = this->mDiff*chebVal.topRows(this->mspSetup->specSize());

      // Compute simple projection
      } else
      {
         // Copy into other array
         this->mTmpZ.topRows(this->mspSetup->specSize()).real() = chebVal.topRows(this->mspSetup->specSize());
      }

      // Set the padded values to zero
      this->mTmpZ.imag().setZero();
      this->mTmpZ.bottomRows(this->mspSetup->padSize()+1).real().setZero();
      this->mTmpZ = this->mPhase.asDiagonal()*this->mTmpZ;

      // Do transform
      checkCudaErrors(cudaMemcpy(this->mDevZ, this->mTmpZ.data(), sizeof(cufftDoubleComplex)*this->mTmpZ.size(), cudaMemcpyHostToDevice));
      checkCudaErrors(cufftExecZ2D(this->mBPlan, this->mDevZ, this->mDevR));
      checkCudaErrors(cudaMemcpy(this->mTmpR.data(), this->mDevR, sizeof(cufftDoubleReal)*this->mTmpR.size(), cudaMemcpyDeviceToHost));

      // Extract physical values
      rPhysVal = this->mTmpR.topRows(this->mspSetup->fwdSize());
   }

   void ChebyshevCuFftTransform::integrate(MatrixZ& rChebVal, const MatrixZ& physVal, ChebyshevCuFftTransform::IntegratorType::Id integrator)
   {
      // Assert that a mixed transform was setup
      assert(this->mspSetup->type() == FftSetup::COMPONENT);

      // assert right sizes for input matrix
      assert(physVal.rows() == this->mspSetup->fwdSize());
      assert(physVal.cols() == this->mspSetup->howmany());

      // assert right sizes for output matrix
      assert(rChebVal.rows() == this->mspSetup->bwdSize());
      assert(rChebVal.cols() == this->mspSetup->howmany());

      // Do transform of real part
      this->mTmpR.topRows(this->mspSetup->fwdSize()) = physVal.real().colwise().reverse();
      this->mTmpR.bottomRows(this->mspSetup->fwdSize()) = physVal.real();
      checkCudaErrors(cudaMemcpy(this->mDevR, this->mTmpR.data(), sizeof(cufftDoubleReal)*this->mTmpR.size(), cudaMemcpyHostToDevice));
      checkCudaErrors(cufftExecD2Z(this->mFPlan, this->mDevR, this->mDevZ));
      checkCudaErrors(cudaMemcpy(this->mTmpZ.data(), this->mDevZ, sizeof(cufftDoubleComplex)*this->mTmpZ.size(), cudaMemcpyDeviceToHost));

      // Rescale FFT output
      this->mTmpZ = this->mPhase_1.asDiagonal()*this->mTmpZ;
      rChebVal.real() = this->mspSetup->scale()*this->mTmpZ.topRows(this->mspSetup->bwdSize()).real();

      // Do transform of imaginary part
      this->mTmpR.topRows(this->mspSetup->fwdSize()) = physVal.imag().colwise().reverse();
      this->mTmpR.bottomRows(this->mspSetup->fwdSize()) = physVal.imag();
      checkCudaErrors(cudaMemcpy(this->mDevR, this->mTmpR.data(), sizeof(cufftDoubleReal)*this->mTmpR.size(), cudaMemcpyHostToDevice));
      checkCudaErrors(cufftExecD2Z(this->mFPlan, this->mDevR, this->mDevZ));
      checkCudaErrors(cudaMemcpy(this->mTmpZ.data(), this->mDevZ, sizeof(cufftDoubleComplex)*this->mTmpZ.size(), cudaMemcpyDeviceToHost));

      // Rescale FFT output
      this->mTmpZ = this->mPhase_1.asDiagonal()*this->mTmpZ;
      rChebVal.imag() = this->mspSetup->scale()*this->mTmpZ.topRows(this->mspSetup->bwdSize()).real();
   }

   void ChebyshevCuFftTransform::project(MatrixZ& rPhysVal, const MatrixZ& chebVal, ChebyshevCuFftTransform::ProjectorType::Id projector)
   {
      // Assert that a mixed transform was setup
      assert(this->mspSetup->type() == FftSetup::COMPONENT);

      // assert on the padding size
      assert(this->mspSetup->padSize() >= 0);
      assert(this->mspSetup->bwdSize() - this->mspSetup->padSize() >= 0);

      // assert right sizes for input  matrix
      assert(chebVal.rows() == this->mspSetup->bwdSize());
      assert(chebVal.cols() == this->mspSetup->howmany());

      // assert right sizes for output matrix
      assert(rPhysVal.rows() == this->mspSetup->fwdSize());
      assert(rPhysVal.cols() == this->mspSetup->howmany());

      // Compute first derivative of real part
      if(projector == ChebyshevCuFftTransform::ProjectorType::DIFF)
      {
         this->mTmpZ.topRows(this->mspSetup->specSize()).real() = this->mDiff*chebVal.topRows(this->mspSetup->specSize()).real();

      // Compute simple projection of real part
      } else
      {
         // Copy values into simple matrix
         this->mTmpZ.topRows(this->mspSetup->specSize()).real() = chebVal.topRows(this->mspSetup->specSize()).real();
      }

      // Set the padded values to zero
      this->mTmpZ.imag().setZero();
      this->mTmpZ.bottomRows(this->mspSetup->padSize()+1).real().setZero();
      this->mTmpZ = this->mPhase.asDiagonal()*this->mTmpZ;

      // Do transform of real part
      checkCudaErrors(cudaMemcpy(this->mDevZ, this->mTmpZ.data(), sizeof(cufftDoubleComplex)*this->mTmpZ.size(), cudaMemcpyHostToDevice));
      checkCudaErrors(cufftExecZ2D(this->mBPlan, this->mDevZ, this->mDevR));
      checkCudaErrors(cudaMemcpy(this->mTmpR.data(), this->mDevR, sizeof(cufftDoubleReal)*this->mTmpR.size(), cudaMemcpyDeviceToHost));
      rPhysVal.real() = this->mTmpR.bottomRows(this->mspSetup->fwdSize());

      // Compute first derivative of imaginary part
      if(projector == ChebyshevCuFftTransform::ProjectorType::DIFF)
      {
         this->mTmpZ.topRows(this->mspSetup->specSize()).real() = this->mDiff*chebVal.topRows(this->mspSetup->specSize()).imag();

      // Compute simple projection of imaginary part
      } else
      {
         // Rescale results
         this->mTmpZ.topRows(this->mspSetup->specSize()).real() = chebVal.topRows(this->mspSetup->specSize()).imag();
      }

      // Set the padded values to zero
      this->mTmpZ.imag().setZero();
      this->mTmpZ.bottomRows(this->mspSetup->padSize()+1).real().setZero();
      this->mTmpZ = this->mPhase.asDiagonal()*this->mTmpZ;

      // Do transform of imaginary part
      checkCudaErrors(cudaMemcpy(this->mDevZ, this->mTmpZ.data(), sizeof(cufftDoubleComplex)*this->mTmpZ.size(), cudaMemcpyHostToDevice));
      checkCudaErrors(cufftExecZ2D(this->mBPlan, this->mDevZ, this->mDevR));
      checkCudaErrors(cudaMemcpy(this->mTmpR.data(), this->mDevR, sizeof(cufftDoubleReal)*this->mTmpR.size(), cudaMemcpyDeviceToHost));
      rPhysVal.imag() = this->mTmpR.bottomRows(this->mspSetup->fwdSize());
   }

#ifdef QUICC_STORAGEPROFILE
   MHDFloat ChebyshevCuFftTransform::requiredStorage() const
   {
      MHDFloat mem = 0.0;

      // Storage required for the fftw plans 
      mem += 8.0*2.0;

      return mem;
   }
#endif // QUICC_STORAGEPROFILE

}
}
