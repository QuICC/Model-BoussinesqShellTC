/** 
 * @file CylinderChebyshevFftwTransform.cpp
 * @brief Source of the implementation of the Chebyshev FFTW transform for a cylinder radius
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "FastTransforms/CylinderChebyshevFftwTransform.hpp"

// Project includes
//
#include "Exceptions/Exception.hpp"
#include "Base/MathConstants.hpp"
#include "FastTransforms/FftwLibrary.hpp"
#include "Python/PythonWrapper.hpp"

namespace GeoMHDiSCC {

namespace Transform {

   Array CylinderChebyshevFftwTransform::generateGrid(const int size)
   {
      // Initialise grid storage
      Array grid(size);

      // Create Chebyshev grid
      for(int k = 0; k < size; k++)
      {
         grid(k) = std::cos((Math::PI)*(static_cast<MHDFloat>(k)+0.5)/static_cast<MHDFloat>(2*size));
      }

      return grid;
   }

   CylinderChebyshevFftwTransform::CylinderChebyshevFftwTransform()
      : mFEPlan(NULL), mFBOPlan(NULL), mBEPlan(NULL)

   {
      PythonWrapper::init();
   }

   CylinderChebyshevFftwTransform::~CylinderChebyshevFftwTransform()
   {
      // Cleanup memory used by FFTW
      FftwLibrary::cleanupFft();
   }

   void CylinderChebyshevFftwTransform::init(CylinderChebyshevFftwTransform::SharedSetupType spSetup)
   {
      // Store the shared pointer to setup object
      this->mspSetup = spSetup;

      // Set the scaling factor
      this->mspSetup->setScale(1.0/static_cast<MHDFloat>(2*this->mspSetup->fwdSize()));

      // Initialise FFTW interface
      this->initFft();

      // Initialise Chebyshev operator(s)
      this->initOperators();

      // Register the FFTW object
      FftwLibrary::registerFft();
   }

   void CylinderChebyshevFftwTransform::requiredOptions(std::set<NonDimensional::Id>& list, const Dimensions::Transform::Id dimId) const
   {
      //
      // No possible options
      //
   }

   void CylinderChebyshevFftwTransform::setOptions(const std::map<NonDimensional::Id, MHDFloat>& options, const Dimensions::Transform::Id dimId)
   {
      //
      // No possible options
      //
   }

   Array CylinderChebyshevFftwTransform::meshGrid() const
   {
      return CylinderChebyshevFftwTransform::generateGrid(this->mspSetup->fwdSize());
   }

   void CylinderChebyshevFftwTransform::initFft()
   {  
      /// \mhdBug implement strided tranforms for complex <-> complex case if possible

      int fwdSize = this->mspSetup->fwdSize();
      int bwdSize = this->mspSetup->bwdSize();
      int howmanyEven = this->mspSetup->howmany(0);
      int howmanyOdd = this->mspSetup->howmany(1);

      // Create the two plans
      const int  *fftSize = &fwdSize;

      if(this->mspSetup->type() == FftSetup::COMPONENT)
      {
      } else
      {
      }

      // Safety assert
      assert(fwdSize == bwdSize);

      // Initialise temporary storage
      this->mTmpEIn.setZero(fwdSize, howmanyEven);
      this->mTmpEOut.setZero(bwdSize, howmanyEven);

      // Create the physical to spectral plan
      const fftw_r2r_kind fwdKind[] = {FFTW_REDFT10};
      this->mFEPlan = fftw_plan_many_r2r(1, fftSize, howmanyEven, this->mTmpEIn.data(), NULL, 1, fwdSize, this->mTmpEOut.data(), NULL, 1, bwdSize, fwdKind, FftwLibrary::planFlag());

      // Create the spectral to physical plan
      const fftw_r2r_kind bwdKind[] = {FFTW_REDFT01};
      this->mBEPlan = fftw_plan_many_r2r(1, fftSize, howmanyEven, this->mTmpEOut.data(), NULL, 1, bwdSize, this->mTmpEIn.data(), NULL, 1, fwdSize, bwdKind, FftwLibrary::planFlag());

      // Initialise temporary storage
      this->mTmpOIn.setZero(fwdSize, howmanyOdd);
      this->mTmpOOut.setZero(bwdSize, howmanyOdd);

      const fftw_r2r_kind fbwdKind[] = {FFTW_REDFT11};
      this->mFBOPlan = fftw_plan_many_r2r(1, fftSize, howmanyOdd, this->mTmpOIn.data(), NULL, 1, fwdSize, this->mTmpOOut.data(), NULL, 1, bwdSize, fbwdKind, FftwLibrary::planFlag());

      if(howmanyEven != howmanyOdd)
      {
         // Initialise temporary storage
         this->mTmpOIn.setZero(fwdSize, howmanyOdd);
         this->mTmpOOut.setZero(bwdSize, howmanyOdd);

         // Create the physical to spectral plan
         const fftw_r2r_kind fwdKind[] = {FFTW_REDFT10};
         this->mFEOPlan = fftw_plan_many_r2r(1, fftSize, howmanyOdd, this->mTmpOIn.data(), NULL, 1, fwdSize, this->mTmpOOut.data(), NULL, 1, bwdSize, fwdKind, FftwLibrary::planFlag());

         // Initialise temporary storage
         this->mTmpOIn.setZero(fwdSize, howmanyOdd);
         this->mTmpOOut.setZero(bwdSize, howmanyOdd);

         // Create the spectral to physical plan
         const fftw_r2r_kind bwdKind[] = {FFTW_REDFT01};
         this->mBEOPlan = fftw_plan_many_r2r(1, fftSize, howmanyOdd, this->mTmpOOut.data(), NULL, 1, bwdSize, this->mTmpOIn.data(), NULL, 1, fwdSize, bwdKind, FftwLibrary::planFlag());

         // Initialise temporary storage
         this->mTmpEIn.setZero(fwdSize, howmanyEven);
         this->mTmpEOut.setZero(bwdSize, howmanyEven);

         const fftw_r2r_kind fbwdKind[] = {FFTW_REDFT11};
         this->mFBOEPlan = fftw_plan_many_r2r(1, fftSize, howmanyEven, this->mTmpOIn.data(), NULL, 1, fwdSize, this->mTmpOOut.data(), NULL, 1, bwdSize, fbwdKind, FftwLibrary::planFlag());
      } else
      {
         this->mFEOPlan = this->mFEPlan;

         this->mBEOPlan = this->mBEPlan;

         this->mFBOEPlan = this->mFBOPlan;
      }
   }

   void CylinderChebyshevFftwTransform::initOperators()
   {
      // First derivative
      this->mDiffE.resize(this->mspSetup->specSize(),this->mspSetup->specSize());
      this->mDiffO.resize(this->mspSetup->specSize(),this->mspSetup->specSize());

      // Initialise python wrapper
      PythonWrapper::import("geomhdiscc.geometry.cylindrical.cylinder_radius");

      #if defined GEOMHDISCC_TRANSOP_FORWARD
         // Initialise array for division by R
         this->mDivR = this->meshGrid().array().pow(-1);

         // Initialise array for division by R^2
         this->mDivR2 = this->meshGrid().array().pow(-2);

         // Prepare arguments to d1(...) call
         PyObject *pArgs, *pValue;
         pArgs = PyTuple_New(3);
         // ... get operator size
         pValue = PyLong_FromLong(this->mspSetup->specSize());
         PyTuple_SetItem(pArgs, 0, pValue);
         // ... create boundray condition (none)
         pValue = PyDict_New();
         PyDict_SetItem(pValue, PyLong_FromLong(0), PyLong_FromLong(0));
         PyTuple_SetItem(pArgs, 2, pValue);

         // ... set even parity
         pValue = PyLong_FromLong(0);
         PyTuple_SetItem(pArgs, 1, pValue);

         // Call d1
         PythonWrapper::setFunction("d1");
         pValue = PythonWrapper::callFunction(pArgs);

         // Fill matrix and clenup
         PythonWrapper::fillMatrix(this->mDiffE, pValue);
         Py_DECREF(pValue);

         // ... set odd parity
         pValue = PyLong_FromLong(1);
         PyTuple_SetItem(pArgs, 1, pValue);

         // Call d1
         PythonWrapper::setFunction("d1");
         pValue = PythonWrapper::callFunction(pArgs);

         // Fill matrix
         PythonWrapper::fillMatrix(this->mDiffO, pValue);
         Py_DECREF(pValue);
      #elif defined GEOMHDISCC_TRANSOP_BACKWARD
         // Division by R
         this->mDivRE.resize(this->mspSetup->specSize(),this->mspSetup->specSize());
         this->mDivRO.resize(this->mspSetup->specSize(),this->mspSetup->specSize()); 

         // Division by R2
         this->mDivR2E.resize(this->mspSetup->specSize(),this->mspSetup->specSize());
         this->mDivR2O.resize(this->mspSetup->specSize(),this->mspSetup->specSize()); 

         // Prepare arguments to i1(...) call
         PyObject *pArgs, *pValue;
         pArgs = PyTuple_New(3);
         // ... get operator size
         pValue = PyLong_FromLong(this->mspSetup->specSize());
         PyTuple_SetItem(pArgs, 0, pValue);
         // ... create boundray condition (last mode zero)
         pValue = PyDict_New();
         PyDict_SetItem(pValue, PyLong_FromLong(0), PyLong_FromLong(99));
         PyTuple_SetItem(pArgs, 2, pValue);

         // ... set even parity
         pValue = PyLong_FromLong(0);
         PyTuple_SetItem(pArgs, 1, pValue);

         // Call i1
         PythonWrapper::setFunction("i1");
         pValue = PythonWrapper::callFunction(pArgs);

         // Fill matrix and clenup
         PythonWrapper::fillMatrix(this->mDiffE, pValue);
         Py_DECREF(pValue);

         // ... set odd parity
         pValue = PyLong_FromLong(1);
         PyTuple_SetItem(pArgs, 1, pValue);

         // Call i1
         PythonWrapper::setFunction("i1");
         pValue = PythonWrapper::callFunction(pArgs);

         // Fill matrix
         PythonWrapper::fillMatrix(this->mDiffO, pValue);
         Py_DECREF(pValue);

         // ... set even parity
         pValue = PyLong_FromLong(0);
         PyTuple_SetItem(pArgs, 1, pValue);

         // Call x1
         PythonWrapper::setFunction("x1");
         pValue = PythonWrapper::callFunction(pArgs);

         // Fill matrix and clenup
         PythonWrapper::fillMatrix(this->mDivRE, pValue);
         Py_DECREF(pValue);

         // ... set odd parity
         pValue = PyLong_FromLong(1);
         PyTuple_SetItem(pArgs, 1, pValue);

         // Call x1
         PythonWrapper::setFunction("x1");
         pValue = PythonWrapper::callFunction(pArgs);

         // Fill matrix
         PythonWrapper::fillMatrix(this->mDivRO, pValue);
         Py_DECREF(pValue);

         // ... set even parity
         pValue = PyLong_FromLong(0);
         PyTuple_SetItem(pArgs, 1, pValue);

         // Call x2
         PythonWrapper::setFunction("x2");
         pValue = PythonWrapper::callFunction(pArgs);

         // Fill matrix and clenup
         PythonWrapper::fillMatrix(this->mDivR2E, pValue);
         Py_DECREF(pValue);

         // ... set odd parity
         pValue = PyLong_FromLong(1);
         PyTuple_SetItem(pArgs, 1, pValue);

         // Call x2
         PythonWrapper::setFunction("x2");
         pValue = PythonWrapper::callFunction(pArgs);

         // Fill matrix
         PythonWrapper::fillMatrix(this->mDivR2O, pValue);
         Py_DECREF(pValue);
      #endif //defined GEOMHDISCC_TRANSOP_FORWARD

      // Fill matrix and cleanup
      PythonWrapper::finalize();

      #if defined GEOMHDISCC_TRANSOP_BACKWARD
         // Factorize matrix and free memory
         this->mSDiffE.compute(this->mDiffE);
         // Check for successful factorisation
         if(this->mSDiffE.info() != Eigen::Success)
         {
            throw Exception("Factorization of backward even differentiation failed!");
         }
         // Factorize matrix and free memory
         this->mSDiffO.compute(this->mDiffO);
         // Check for successful factorisation
         if(this->mSDiffO.info() != Eigen::Success)
         {
            throw Exception("Factorization of backward odd differentiation failed!");
         }

         // Factorize division matrix and free memory
         this->mSDivRE.compute(this->mDivRE);
         // Check for successful factorisation
         if(this->mSDivRE.info() != Eigen::Success)
         {
            throw Exception("Factorization of backward even division by R failed!");
         }
         // Factorize division matrix and free memory
         this->mSDivRO.compute(this->mDivRO);
         // Check for successful factorisation
         if(this->mSDivRO.info() != Eigen::Success)
         {
            throw Exception("Factorization of backward odd division by R failed!");
         }

         // Factorize division matrix and free memory
         this->mSDivR2E.compute(this->mDivR2E);
         // Check for successful factorisation
         if(this->mSDivR2E.info() != Eigen::Success)
         {
            throw Exception("Factorization of backward even division by R^2 failed!");
         }
         // Factorize division matrix and free memory
         this->mSDivR2O.compute(this->mDivR2O);
         // Check for successful factorisation
         if(this->mSDivR2O.info() != Eigen::Success)
         {
            throw Exception("Factorization of backward odd division by R^2 failed!");
         }
      #endif //defined GEOMHDISCC_TRANSOP_BACKWARD
   }

   void CylinderChebyshevFftwTransform::cleanupFft()
   {
      // Detroy forward even plan
      if(this->mFEPlan)
      {
         fftw_destroy_plan(this->mFEPlan);
      }

      // Detroy backward even plan
      if(this->mBEPlan)
      {
         fftw_destroy_plan(this->mBEPlan);
      }

      // Detroy forward/backward odd plan
      if(this->mFBOPlan)
      {
         fftw_destroy_plan(this->mFBOPlan);
      }

      // Detroy forward even plan on odd sizes
      if(this->mFEOPlan)
      {
         fftw_destroy_plan(this->mFEOPlan);
      }

      // Detroy backward even plan on odd sizes
      if(this->mBEOPlan)
      {
         fftw_destroy_plan(this->mBEOPlan);
      }

      // Detroy forward/backward odd plan on even sizes
      if(this->mFBOEPlan)
      {
         fftw_destroy_plan(this->mFBOEPlan);
      }

      // Unregister the FFTW object
      FftwLibrary::unregisterFft();

      // cleanup fftw library
      FftwLibrary::cleanupFft();
   }

   void CylinderChebyshevFftwTransform::extractParityModes(Matrix& rSelected, const Matrix& data, const MatrixI& info, const int rows)
   {
      assert(rSelected.cols() == info.col(1).sum());

      int k = 0;
      for(int i = 0; i < info.rows(); ++i)
      {
         int j0 = info(i,0);
         for(int j = 0; j < info(i,1); ++j, ++k)
         { 
            rSelected.col(k).topRows(rows) = data.col(j0 + j).topRows(rows);
         }
      }
   }

   void CylinderChebyshevFftwTransform::setParityModes(Matrix& rData, const Matrix& selected, const MatrixI& info, const int rows)
   {
      assert(selected.cols() == info.col(1).sum());

      int k = 0;
      for(int i = 0; i < info.rows(); ++i)
      {
         int j0 = info(i,0);
         for(int j = 0; j < info(i,1); ++j, ++k)
         { 
            rData.col(j0 + j).topRows(rows) = selected.col(k).topRows(rows);
         }
      }
   }

   void CylinderChebyshevFftwTransform::addParityModes(Matrix& rData, const Matrix& selected, const MatrixI& info, const int rows)
   {
      assert(selected.cols() == info.col(1).sum());

      int k = 0;
      for(int i = 0; i < info.rows(); ++i)
      {
         int j0 = info(i,0);
         for(int j = 0; j < info(i,1); ++j, ++k)
         { 
            rData.col(j0 + j).topRows(rows) += selected.col(k).topRows(rows);
         }
      }
   }

   void CylinderChebyshevFftwTransform::extractParityModes(Matrix& rSelected, const MatrixZ& data, const bool isReal, const MatrixI& info, const int rows)
   {
      assert(rSelected.cols() == info.col(1).sum());

      int k = 0;
      if(isReal)
      {
         for(int i = 0; i < info.rows(); ++i)
         {
            int j0 = info(i,0);
            for(int j = 0; j < info(i,1); ++j, ++k)
            { 
               rSelected.col(k).topRows(rows) = data.col(j0 + j).topRows(rows).real();
            }
         }
      } else
      {
         for(int i = 0; i < info.rows(); ++i)
         {
            int j0 = info(i,0);
            for(int j = 0; j < info(i,1); ++j, ++k)
            { 
               rSelected.col(k).topRows(rows) = data.col(j0 + j).topRows(rows).imag();
            }
         }
      }
   }

   void CylinderChebyshevFftwTransform::setParityModes(MatrixZ& rData, const Matrix& selected, const bool isReal, const MatrixI& info, const int rows)
   {
      assert(selected.cols() == info.col(1).sum());

      int k = 0;
      if(isReal)
      {
         for(int i = 0; i < info.rows(); ++i)
         {
            int j0 = info(i,0);
            for(int j = 0; j < info(i,1); ++j, ++k)
            { 
               rData.col(j0 + j).topRows(rows).real() = selected.col(k).topRows(rows);
            }
         }
      } else
      {
         for(int i = 0; i < info.rows(); ++i)
         {
            int j0 = info(i,0);
            for(int j = 0; j < info(i,1); ++j, ++k)
            { 
               rData.col(j0 + j).topRows(rows).imag() = selected.col(k).topRows(rows);
            }
         }
      }
   }

   void CylinderChebyshevFftwTransform::addParityModes(MatrixZ& rData, const Matrix& selected, const bool isReal, const MatrixI& info, const int rows)
   {
      assert(selected.cols() == info.col(1).sum());

      int k = 0;
      if(isReal)
      {
         for(int i = 0; i < info.rows(); ++i)
         {
            int j0 = info(i,0);
            for(int j = 0; j < info(i,1); ++j, ++k)
            { 
               rData.col(j0 + j).topRows(rows).real() += selected.col(k).topRows(rows);
            }
         }
      } else
      {
         for(int i = 0; i < info.rows(); ++i)
         {
            int j0 = info(i,0);
            for(int j = 0; j < info(i,1); ++j, ++k)
            { 
               rData.col(j0 + j).topRows(rows).imag() += selected.col(k).topRows(rows);
            }
         }
      }
   }

#ifdef GEOMHDISCC_STORAGEPROFILE
   MHDFloat CylinderChebyshevFftwTransform::requiredStorage() const
   {
      MHDFloat mem = 0.0;

      // Storage required for the fftw plans 
      mem += 8.0*2.0;

      return mem;
   }
#endif // GEOMHDISCC_STORAGEPROFILE

}
}
