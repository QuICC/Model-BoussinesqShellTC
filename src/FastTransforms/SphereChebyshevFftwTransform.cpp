/** 
 * @file SphereChebyshevFftwTransform.cpp
 * @brief Source of the implementation of the Chebyshev FFTW transform for a sphere radius
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// System includes
//
#include <cassert>
#include <limits>

// External includes
//

// Class include
//
#include "FastTransforms/SphereChebyshevFftwTransform.hpp"

// Project includes
//
#include "StaticAsserts/StaticAssert.hpp"
#include "Exceptions/Exception.hpp"
#include "Base/MathConstants.hpp"
#include "FastTransforms/FftwLibrary.hpp"
#include "FastTransforms/ParityTransformTools.hpp"
#include "Python/PythonWrapper.hpp"

#include <iostream>
namespace GeoMHDiSCC {

namespace Transform {

   Array SphereChebyshevFftwTransform::generateGrid(const int size)
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

   SphereChebyshevFftwTransform::SphereChebyshevFftwTransform()
      : mFEPlan(NULL), mFBOPlan(NULL), mBEPlan(NULL)
   {
      // Initialise the Python interpreter wrapper
      PythonWrapper::init();

      // Initialize FFTW
      FftwLibrary::initFft();
   }

   SphereChebyshevFftwTransform::~SphereChebyshevFftwTransform()
   {
      // Cleanup memory used by FFTW
      FftwLibrary::cleanupFft();
   }

   void SphereChebyshevFftwTransform::init(SphereChebyshevFftwTransform::SharedSetupType spSetup)
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

   void SphereChebyshevFftwTransform::requiredOptions(std::set<NonDimensional::Id>& list, const Dimensions::Transform::Id dimId) const
   {
      //
      // No possible options
      //
   }

   void SphereChebyshevFftwTransform::setOptions(const std::map<NonDimensional::Id, MHDFloat>& options, const Dimensions::Transform::Id dimId)
   {
      //
      // No possible options
      //
   }

   Array SphereChebyshevFftwTransform::meshGrid() const
   {
      return SphereChebyshevFftwTransform::generateGrid(this->mspSetup->fwdSize());
   }

   void SphereChebyshevFftwTransform::initFft()
   {  
      /// \mhdBug implement strided tranforms for complex <-> complex case if possible

      int fwdSize = this->mspSetup->fwdSize();
      int bwdSize = this->mspSetup->bwdSize();
      int howmanyEven = this->mspSetup->howmany(0);
      int howmanyOdd = this->mspSetup->howmany(1);
      int howmany = std::max(howmanyEven, howmanyOdd);

      // Create the two plans
      const int  *fftSize = &fwdSize;

      if(this->mspSetup->type() == FftSetup::COMPONENT)
      {
      } else
      {
      }

      // Safety assert
      assert(fwdSize == bwdSize);

      // Initialise temporary storage (to max size)
      this->mTmpIn.setZero(fwdSize, howmany);
      this->mTmpOut.setZero(bwdSize, howmany);

      // Create the physical to spectral plan
      const fftw_r2r_kind fwdKind[] = {FFTW_REDFT10};
      this->mFEPlan = fftw_plan_many_r2r(1, fftSize, howmanyEven, this->mTmpIn.data(), NULL, 1, fwdSize, this->mTmpOut.data(), NULL, 1, bwdSize, fwdKind, FftwLibrary::planFlag());

      // Create the spectral to physical plan
      const fftw_r2r_kind bwdKind[] = {FFTW_REDFT01};
      this->mBEPlan = fftw_plan_many_r2r(1, fftSize, howmanyEven, this->mTmpOut.data(), NULL, 1, bwdSize, this->mTmpIn.data(), NULL, 1, fwdSize, bwdKind, FftwLibrary::planFlag());

      const fftw_r2r_kind fbwdKind[] = {FFTW_REDFT11};
      this->mFBOPlan = fftw_plan_many_r2r(1, fftSize, howmanyOdd, this->mTmpIn.data(), NULL, 1, fwdSize, this->mTmpOut.data(), NULL, 1, bwdSize, fbwdKind, FftwLibrary::planFlag());
      
      // Different plans are needed when parity changes!
      if(howmanyEven != howmanyOdd)
      {
         // Create the physical to spectral plan
         const fftw_r2r_kind fwdKind[] = {FFTW_REDFT10};
         this->mFEOPlan = fftw_plan_many_r2r(1, fftSize, howmanyOdd, this->mTmpIn.data(), NULL, 1, fwdSize, this->mTmpOut.data(), NULL, 1, bwdSize, fwdKind, FftwLibrary::planFlag());

         // Create the spectral to physical plan
         const fftw_r2r_kind bwdKind[] = {FFTW_REDFT01};
         this->mBEOPlan = fftw_plan_many_r2r(1, fftSize, howmanyOdd, this->mTmpOut.data(), NULL, 1, bwdSize, this->mTmpIn.data(), NULL, 1, fwdSize, bwdKind, FftwLibrary::planFlag());

         const fftw_r2r_kind fbwdKind[] = {FFTW_REDFT11};
         this->mFBOEPlan = fftw_plan_many_r2r(1, fftSize, howmanyEven, this->mTmpIn.data(), NULL, 1, fwdSize, this->mTmpOut.data(), NULL, 1, bwdSize, fbwdKind, FftwLibrary::planFlag());
      } else
      {
         this->mFEOPlan = this->mFEPlan;

         this->mBEOPlan = this->mBEPlan;

         this->mFBOEPlan = this->mFBOPlan;
      }
   }

   void SphereChebyshevFftwTransform::initOperators()
   {
      //
      // Initialise projector operators
      //
      this->mProjectorFlips.insert(std::make_pair(ProjectorType::PROJ, 0));
      this->mProjectorFlips.insert(std::make_pair(ProjectorType::DIVR, 1));
      this->mProjectorFlips.insert(std::make_pair(ProjectorType::DIVR2, 0));
      this->mProjectorFlips.insert(std::make_pair(ProjectorType::DIFF, 1));
      this->mProjectorFlips.insert(std::make_pair(ProjectorType::DIFF2, 0));
      this->mProjectorFlips.insert(std::make_pair(ProjectorType::DIFFR, 0));
      this->mProjectorFlips.insert(std::make_pair(ProjectorType::DIVRDIFFR, 1));
      this->mProjectorFlips.insert(std::make_pair(ProjectorType::RADLAPL, 0));

      //
      // Initialise integrator operators
      //
      this->mIntegratorFlips.insert(std::make_pair(IntegratorType::INTG, 0));
      std::pair<SparseMatrix,SparseMatrix> opPair = std::make_pair(SparseMatrix(this->mspSetup->fwdSize(),this->mspSetup->fwdSize()),SparseMatrix(this->mspSetup->fwdSize(),this->mspSetup->fwdSize()));
      // Multiplication by R
      this->mIntgOp.insert(std::make_pair(IntegratorType::INTGR, opPair));
      this->mIntegratorFlips.insert(std::make_pair(IntegratorType::INTGR, 1));
      // QST Q operator (4th order)
      this->mIntgOp.insert(std::make_pair(IntegratorType::INTGQ4, opPair));
      this->mIntegratorFlips.insert(std::make_pair(IntegratorType::INTGQ4, 1));
      // QST S operator (4th order)
      this->mIntgOp.insert(std::make_pair(IntegratorType::INTGS4, opPair));
      this->mIntegratorFlips.insert(std::make_pair(IntegratorType::INTGS4, 1));
      // QST T operator
      this->mIntgOp.insert(std::make_pair(IntegratorType::INTGT, opPair));
      this->mIntegratorFlips.insert(std::make_pair(IntegratorType::INTGT, 0));
      // QST Q operator (2th order)
      this->mIntgOp.insert(std::make_pair(IntegratorType::INTGQ2, opPair));
      this->mIntegratorFlips.insert(std::make_pair(IntegratorType::INTGQ2, 1));
      // QST S operator (2th order)
      this->mIntgOp.insert(std::make_pair(IntegratorType::INTGS2, opPair));
      this->mIntegratorFlips.insert(std::make_pair(IntegratorType::INTGS2, 1));

      //
      // Initialise regularity operators
      //
      // 0th order Regularity constraint
      this->mRegOp.insert(std::make_pair(RegularityType::REG0, opPair));
      this->mRegularitySize.insert(std::make_pair(RegularityType::REG0, std::make_pair(1,0)));
      // 1st order Regularity constraint
      this->mRegOp.insert(std::make_pair(RegularityType::REG1, opPair));
      this->mRegularitySize.insert(std::make_pair(RegularityType::REG1, std::make_pair(1,1)));
      // 2nd order Regularity constraint
      this->mRegOp.insert(std::make_pair(RegularityType::REG2, opPair));
      this->mRegularitySize.insert(std::make_pair(RegularityType::REG2, std::make_pair(2,1)));
      // 3rd order Regularity constraint
      this->mRegOp.insert(std::make_pair(RegularityType::REG3, opPair));
      this->mRegularitySize.insert(std::make_pair(RegularityType::REG3, std::make_pair(2,2)));
      // 4th order Regularity constraint
      this->mRegOp.insert(std::make_pair(RegularityType::REG4, opPair));
      this->mRegularitySize.insert(std::make_pair(RegularityType::REG4, std::make_pair(3,2)));

      //
      // Initialise solver operators
      //
      // Multiplication by R
      this->mSolveOp.insert(std::make_pair(ProjectorType::DIVR, opPair));
      // Multiplication by R^2
      this->mSolveOp.insert(std::make_pair(ProjectorType::DIVR2, opPair));
      // First derivative
      this->mSolveOp.insert(std::make_pair(ProjectorType::DIFF, opPair));
      // Second derivative
      this->mSolveOp.insert(std::make_pair(ProjectorType::DIFF2, opPair));

      //
      // Initialise regularity solver operators
      //
      // 0th order Regularity constraint
      this->mRegSolveOp.insert(std::make_pair(RegularityType::REG0, opPair));
      // 1st order Regularity constraint
      this->mRegSolveOp.insert(std::make_pair(RegularityType::REG1, opPair));
      // 2nd order Regularity constraint
      this->mRegSolveOp.insert(std::make_pair(RegularityType::REG2, opPair));
      // 3rd order Regularity constraint
      this->mRegSolveOp.insert(std::make_pair(RegularityType::REG3, opPair));
      // 4th order Regularity constraint
      this->mRegSolveOp.insert(std::make_pair(RegularityType::REG4, opPair));

      //
      // Initialise solvers
      //
      SharedPtrMacro<Solver::SparseSelector<SparseMatrix>::Type>  pSolE(new Solver::SparseSelector<SparseMatrix>::Type());
      SharedPtrMacro<Solver::SparseSelector<SparseMatrix>::Type>  pSolO(new Solver::SparseSelector<SparseMatrix>::Type());
      this->mSolver.insert(std::make_pair(ProjectorType::DIVR, std::make_pair(pSolE,pSolO)));
      pSolE = SharedPtrMacro<Solver::SparseSelector<SparseMatrix>::Type>(new Solver::SparseSelector<SparseMatrix>::Type());
      pSolO = SharedPtrMacro<Solver::SparseSelector<SparseMatrix>::Type>(new Solver::SparseSelector<SparseMatrix>::Type());
      this->mSolver.insert(std::make_pair(ProjectorType::DIVR2, std::make_pair(pSolE,pSolO)));
      pSolE = SharedPtrMacro<Solver::SparseSelector<SparseMatrix>::Type>(new Solver::SparseSelector<SparseMatrix>::Type());
      pSolO = SharedPtrMacro<Solver::SparseSelector<SparseMatrix>::Type>(new Solver::SparseSelector<SparseMatrix>::Type());
      this->mSolver.insert(std::make_pair(ProjectorType::DIFF, std::make_pair(pSolE,pSolO)));
      pSolE = SharedPtrMacro<Solver::SparseSelector<SparseMatrix>::Type>(new Solver::SparseSelector<SparseMatrix>::Type());
      pSolO = SharedPtrMacro<Solver::SparseSelector<SparseMatrix>::Type>(new Solver::SparseSelector<SparseMatrix>::Type());
      this->mSolver.insert(std::make_pair(ProjectorType::DIFF2, std::make_pair(pSolE,pSolO)));

      //
      // Initialise regularity solvers
      //
      pSolE = SharedPtrMacro<Solver::SparseSelector<SparseMatrix>::Type>(new Solver::SparseSelector<SparseMatrix>::Type());
      pSolO = SharedPtrMacro<Solver::SparseSelector<SparseMatrix>::Type>(new Solver::SparseSelector<SparseMatrix>::Type());
      this->mRegSolver.insert(std::make_pair(RegularityType::REG0, std::make_pair(pSolE,pSolO)));
      pSolE = SharedPtrMacro<Solver::SparseSelector<SparseMatrix>::Type>(new Solver::SparseSelector<SparseMatrix>::Type());
      pSolO = SharedPtrMacro<Solver::SparseSelector<SparseMatrix>::Type>(new Solver::SparseSelector<SparseMatrix>::Type());
      this->mRegSolver.insert(std::make_pair(RegularityType::REG1, std::make_pair(pSolE,pSolO)));
      pSolE = SharedPtrMacro<Solver::SparseSelector<SparseMatrix>::Type>(new Solver::SparseSelector<SparseMatrix>::Type());
      pSolO = SharedPtrMacro<Solver::SparseSelector<SparseMatrix>::Type>(new Solver::SparseSelector<SparseMatrix>::Type());
      this->mRegSolver.insert(std::make_pair(RegularityType::REG2, std::make_pair(pSolE,pSolO)));
      pSolE = SharedPtrMacro<Solver::SparseSelector<SparseMatrix>::Type>(new Solver::SparseSelector<SparseMatrix>::Type());
      pSolO = SharedPtrMacro<Solver::SparseSelector<SparseMatrix>::Type>(new Solver::SparseSelector<SparseMatrix>::Type());
      this->mRegSolver.insert(std::make_pair(RegularityType::REG3, std::make_pair(pSolE,pSolO)));
      pSolE = SharedPtrMacro<Solver::SparseSelector<SparseMatrix>::Type>(new Solver::SparseSelector<SparseMatrix>::Type());
      pSolO = SharedPtrMacro<Solver::SparseSelector<SparseMatrix>::Type>(new Solver::SparseSelector<SparseMatrix>::Type());
      this->mRegSolver.insert(std::make_pair(RegularityType::REG4, std::make_pair(pSolE,pSolO)));

      //
      // Initialise solvers
      //

      // Initialise python wrapper
      PythonWrapper::import("geomhdiscc.geometry.spherical.sphere_radius");

      // Prepare arguments to Chebyshev matrices call
      PyObject *pArgs, *pValue, *pRegArgs;
      pArgs = PyTuple_New(3);
      pRegArgs = PyTuple_New(4);
      // ... get operator size
      pValue = PyLong_FromLong(this->mspSetup->specSize());
      PyTuple_SetItem(pArgs, 0, pValue);

      for(int parity = 0; parity < 2; ++parity)
      {
         // ... create boundray condition (none)
         pValue = PyDict_New();
         PyDict_SetItem(pValue, PyLong_FromLong(0), PyLong_FromLong(0));
         PyTuple_SetItem(pArgs, 2, pValue);

         // ... set forward operator size
         pValue = PyLong_FromLong(this->mspSetup->fwdSize());
         PyTuple_SetItem(pArgs, 0, pValue);

         // Set parity
         PyTuple_SetItem(pArgs, 1, PyLong_FromLong(parity));

         // Call r1
         PythonWrapper::setFunction("r1");
         pValue = PythonWrapper::callFunction(pArgs);
         // Fill matrix
         PythonWrapper::fillMatrix(this->intgOp(IntegratorType::INTGR, parity), pValue);
         Py_DECREF(pValue);

         // Call QST Q component (4th order)
         PythonWrapper::setFunction("i4r3");
         pValue = PythonWrapper::callFunction(pArgs);
         // Fill matrix
         PythonWrapper::fillMatrix(this->intgOp(IntegratorType::INTGQ4, parity), pValue);
         Py_DECREF(pValue);

         // Call QST S component (4th order)
         PythonWrapper::setFunction("i4r3d1r1");
         pValue = PythonWrapper::callFunction(pArgs);
         // Fill matrix
         PythonWrapper::fillMatrix(this->intgOp(IntegratorType::INTGS4, parity), pValue);
         Py_DECREF(pValue);

         // Call QST T component
         PythonWrapper::setFunction("i2r2");
         pValue = PythonWrapper::callFunction(pArgs);
         // Fill matrix
         PythonWrapper::fillMatrix(this->intgOp(IntegratorType::INTGT, parity), pValue);
         Py_DECREF(pValue);

         // Call QST Q component (2nd order)
         PythonWrapper::setFunction("i2r1");
         pValue = PythonWrapper::callFunction(pArgs);
         // Fill matrix
         PythonWrapper::fillMatrix(this->intgOp(IntegratorType::INTGQ2, parity), pValue);
         Py_DECREF(pValue);

         // Call QST S component (2nd order)
         PythonWrapper::setFunction("i2r1d1r1");
         pValue = PythonWrapper::callFunction(pArgs);
         // Fill matrix
         PythonWrapper::fillMatrix(this->intgOp(IntegratorType::INTGS2, parity), pValue);
         Py_DECREF(pValue);

         // Change resoluton for solver matrices
         pValue = PyLong_FromLong(this->mspSetup->fwdSize());
         PyTuple_SetItem(pArgs, 0, pValue);

         // Call r1 for solver
         PythonWrapper::setFunction("r1");
         pValue = PythonWrapper::callFunction(pArgs);
         // Fill matrix
         PythonWrapper::fillMatrix(this->solveOp(ProjectorType::DIVR, parity), pValue);
         Py_DECREF(pValue);

         // Call r2 for solver
         PythonWrapper::setFunction("r2");
         pValue = PythonWrapper::callFunction(pArgs);
         // Fill matrix
         PythonWrapper::fillMatrix(this->solveOp(ProjectorType::DIVR2, parity), pValue);
         Py_DECREF(pValue);

         // Call d1 for solver
         // ... change boundary condition to zero last modes
         pValue = PyTuple_GetItem(pArgs, 2);
         PyDict_SetItem(pValue, PyLong_FromLong(0), PyLong_FromLong(991));
         PythonWrapper::setFunction("i1");
         pValue = PythonWrapper::callFunction(pArgs);
         // Fill matrix
         PythonWrapper::fillMatrix(this->solveOp(ProjectorType::DIFF, parity), pValue);
         Py_DECREF(pValue);

         // Call d2 for solver
         // ... change boundary condition to zero last modes
         pValue = PyTuple_GetItem(pArgs, 2);
         PyDict_SetItem(pValue, PyLong_FromLong(0), PyLong_FromLong(992));
         PythonWrapper::setFunction("i2");
         pValue = PythonWrapper::callFunction(pArgs);
         // Fill matrix
         PythonWrapper::fillMatrix(this->solveOp(ProjectorType::DIFF2, parity), pValue);
         Py_DECREF(pValue);

         // Prepare operators
         pValue = PyLong_FromLong(this->mspSetup->specSize());
         PyTuple_SetItem(pRegArgs, 0, pValue);
         PyTuple_SetItem(pRegArgs, 1, PyLong_FromLong(parity));
         pValue = PyDict_New();
         PyDict_SetItem(pValue, PyLong_FromLong(0), PyLong_FromLong(0));
         PyTuple_SetItem(pRegArgs, 2, pValue);
         PyTuple_SetItem(pRegArgs, 3, PyLong_FromLong(0));

         if(parity == 0)
         {
            // Call stencil for regularity solver
            // ... change boundary condition to 0th order regularity
            pValue = PyTuple_GetItem(pRegArgs, 2);
            PyDict_SetItem(pValue, PyLong_FromLong(0), PyLong_FromLong(-24));
            PythonWrapper::setFunction("stencil");
            pValue = PythonWrapper::callFunction(pRegArgs);
            // Fill matrix
            PythonWrapper::fillMatrix(this->regOp(RegularityType::REG0, parity), pValue);
            Py_DECREF(pValue);
         }

         // Call stencil for regularity solver
         // ... change boundary condition to 1st order regularity
         pValue = PyTuple_GetItem(pRegArgs, 2);
         PyDict_SetItem(pValue, PyLong_FromLong(0), PyLong_FromLong(-25));
         PythonWrapper::setFunction("stencil");
         pValue = PythonWrapper::callFunction(pRegArgs);
         // Fill matrix
         PythonWrapper::fillMatrix(this->regOp(RegularityType::REG1, parity), pValue);
         Py_DECREF(pValue);

         // Call stencil for regularity solver
         // ... change boundary condition to 2nd order regularity
         pValue = PyTuple_GetItem(pRegArgs, 2);
         PyDict_SetItem(pValue, PyLong_FromLong(0), PyLong_FromLong(-30));
         PythonWrapper::setFunction("stencil");
         pValue = PythonWrapper::callFunction(pRegArgs);
         // Fill matrix
         PythonWrapper::fillMatrix(this->regOp(RegularityType::REG2, parity), pValue);
         Py_DECREF(pValue);

         // Call stencil for regularity solver
         // ... change boundary condition to 3rd order regularity
         pValue = PyTuple_GetItem(pRegArgs, 2);
         PyDict_SetItem(pValue, PyLong_FromLong(0), PyLong_FromLong(-31));
         PythonWrapper::setFunction("stencil");
         pValue = PythonWrapper::callFunction(pRegArgs);
         // Fill matrix
         PythonWrapper::fillMatrix(this->regOp(RegularityType::REG3, parity), pValue);
         Py_DECREF(pValue);

//         // Call stencil for regularity solver
//         // ... change boundary condition to 4th order regularity
//         pValue = PyTuple_GetItem(pRegArgs, 2);
//         PyDict_SetItem(pValue, PyLong_FromLong(0), PyLong_FromLong(-40));
//         PyDict_SetItem(pValue, PyUnicode_FromString("rt"), PyLong_FromLong(3));
//         PythonWrapper::setFunction("stencil");
//         pValue = PythonWrapper::callFunction(pRegArgs);
//         // Fill matrix
//         PythonWrapper::fillMatrix(this->regOp(RegularityType::REG4, parity), pValue);
//         Py_DECREF(pValue);

         // Prepare regularity solvers
         pValue = PyLong_FromLong(this->mspSetup->fwdSize());
         PyTuple_SetItem(pRegArgs, 0, pValue);
         PyTuple_SetItem(pRegArgs, 3, PyLong_FromLong(1));

         // Call stencil for regularity solver
         // ... change boundary condition to 0th order regularity (only even harmonic)
         if(parity == 0)
         {
            pValue = PyTuple_GetItem(pRegArgs, 2);
            PyDict_SetItem(pValue, PyLong_FromLong(0), PyLong_FromLong(-24));
            PyDict_SetItem(pValue, PyUnicode_FromString("rt"), PyLong_FromLong(1-parity));
            PythonWrapper::setFunction("stencil");
            pValue = PythonWrapper::callFunction(pRegArgs);
            // Fill matrix
            PythonWrapper::fillMatrix(this->regSolveOp(RegularityType::REG0, parity), pValue);
            Py_DECREF(pValue);
         }

         // Call stencil for regularity solver
         // ... change boundary condition to 1st order regularity
         pValue = PyTuple_GetItem(pRegArgs, 2);
         PyDict_SetItem(pValue, PyLong_FromLong(0), PyLong_FromLong(-25));
         PyDict_SetItem(pValue, PyUnicode_FromString("rt"), PyLong_FromLong(1));
         PythonWrapper::setFunction("stencil");
         pValue = PythonWrapper::callFunction(pRegArgs);
         // Fill matrix
         PythonWrapper::fillMatrix(this->regSolveOp(RegularityType::REG1, parity), pValue);
         Py_DECREF(pValue);

         // Call stencil for regularity solver
         // ... change boundary condition to 2nd order regularity
         pValue = PyTuple_GetItem(pRegArgs, 2);
         PyDict_SetItem(pValue, PyLong_FromLong(0), PyLong_FromLong(-30));
         PyDict_SetItem(pValue, PyUnicode_FromString("rt"), PyLong_FromLong(2 - parity));
         PythonWrapper::setFunction("stencil");
         pValue = PythonWrapper::callFunction(pRegArgs);
         // Fill matrix
         PythonWrapper::fillMatrix(this->regSolveOp(RegularityType::REG2, parity), pValue);
         Py_DECREF(pValue);

         // Call stencil for regularity solver
         // ... change boundary condition to 3rd order regularity
         pValue = PyTuple_GetItem(pRegArgs, 2);
         PyDict_SetItem(pValue, PyLong_FromLong(0), PyLong_FromLong(-31));
         PyDict_SetItem(pValue, PyUnicode_FromString("rt"), PyLong_FromLong(2));
         PythonWrapper::setFunction("stencil");
         pValue = PythonWrapper::callFunction(pRegArgs);
         // Fill matrix
         PythonWrapper::fillMatrix(this->regSolveOp(RegularityType::REG3, parity), pValue);
         Py_DECREF(pValue);

//         // Call stencil for regularity solver
//         // ... change boundary condition to 4th order regularity
//         pValue = PyTuple_GetItem(pRegArgs, 2);
//         PyDict_SetItem(pValue, PyLong_FromLong(0), PyLong_FromLong(-40));
//         PyDict_SetItem(pValue, PyUnicode_FromString("rt"), PyLong_FromLong(3 - parity));
//         PythonWrapper::setFunction("stencil");
//         pValue = PythonWrapper::callFunction(pRegArgs);
//         // Fill matrix
//         PythonWrapper::fillMatrix(this->regSolveOp(RegularityType::REG4, parity), pValue);
//         Py_DECREF(pValue);

         // Initialize solver storage
         this->tmpInS(parity).setZero(this->mspSetup->fwdSize(), this->mspSetup->howmany(parity));
         this->tmpOutS(parity).setZero(this->mspSetup->bwdSize(), this->mspSetup->howmany(parity));

         // Initialize solver and factorize division by R operator
         this->solver(ProjectorType::DIVR, parity).compute(this->solveOp(ProjectorType::DIVR, parity));
         // Check for successful factorisation
         if(this->solver(ProjectorType::DIVR, parity).info() != Eigen::Success)
         {
            throw Exception("Factorization of backward division by R failed!");
         }

         // Initialize solver and factorize division by R^2 operator
         this->solver(ProjectorType::DIVR2, parity).compute(this->solveOp(ProjectorType::DIVR2, parity));
         // Check for successful factorisation
         if(this->solver(ProjectorType::DIVR2, parity).info() != Eigen::Success)
         {
            throw Exception("Factorization of backward division by R^2 failed!");
         }

         // Initialize solver and factorize division by d1 operator
         this->solver(ProjectorType::DIFF, parity).compute(this->solveOp(ProjectorType::DIFF, parity));
         // Check for successful factorisation
         if(this->solver(ProjectorType::DIFF, parity).info() != Eigen::Success)
         {
            throw Exception("Factorization of backward 1st derivative failed!");
         }

         // Initialize solver and factorize division by d2 operator
         this->solver(ProjectorType::DIFF2, parity).compute(this->solveOp(ProjectorType::DIFF2, parity));
         // Check for successful factorisation
         if(this->solver(ProjectorType::DIFF2, parity).info() != Eigen::Success)
         {
            throw Exception("Factorization of backward 2nd derivative failed!");
         }

         // Initialize solver and factorize 0th order regularity stencil
         if(parity == 0)
         {
            this->regSolver(RegularityType::REG0, parity).compute(this->regSolveOp(RegularityType::REG0, parity));
            // Check for successful factorisation
            if(this->regSolver(RegularityType::REG0, parity).info() != Eigen::Success)
            {
               throw Exception("Factorization of 0th order regularity stencil failed!");
            }
         }

         // Initialize solver and factorize 1st order regularity stencil
         this->regSolver(RegularityType::REG1, parity).compute(this->regSolveOp(RegularityType::REG1, parity));
         // Check for successful factorisation
         if(this->regSolver(RegularityType::REG1, parity).info() != Eigen::Success)
         {
            throw Exception("Factorization of 1st order regularity stencil failed!");
         }

         // Initialize solver and factorize 2nd order regularity stencil
         this->regSolver(RegularityType::REG2, parity).compute(this->regSolveOp(RegularityType::REG2, parity));
         // Check for successful factorisation
         if(this->regSolver(RegularityType::REG2, parity).info() != Eigen::Success)
         {
            throw Exception("Factorization of 2nd order regularity stencil failed!");
         }

         // Initialize solver and factorize 3rd order regularity stencil
         this->regSolver(RegularityType::REG3, parity).compute(this->regSolveOp(RegularityType::REG3, parity));
         // Check for successful factorisation
         if(this->regSolver(RegularityType::REG3, parity).info() != Eigen::Success)
         {
            throw Exception("Factorization of 3rd order regularity stencil failed!");
         }

//         // Initialize solver and factorize 4th order regularity stencil
//         this->regSolver(RegularityType::REG4, parity).compute(this->regSolveOp(RegularityType::REG4, parity));
//         // Check for successful factorisation
//         if(this->regSolver(RegularityType::REG4, parity).info() != Eigen::Success)
//         {
//            throw Exception("Factorization of 4th order regularity stencil failed!");
//         }
      }

      PythonWrapper::finalize();
   }

   void SphereChebyshevFftwTransform::cleanupFft()
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

   void SphereChebyshevFftwTransform::integrate(Matrix& rChebVal, const Matrix& physVal, SphereChebyshevFftwTransform::IntegratorType::Id integrator, Arithmetics::Id arithId)
   {
      assert(arithId == Arithmetics::SET);

      // Assert that a mixed transform was not setup
      assert(this->mspSetup->type() == FftSetup::REAL);

      // assert right sizes for input matrix
      assert(physVal.rows() == this->mspSetup->fwdSize());
      assert(physVal.cols() == this->mspSetup->howmany());

      // assert right sizes for output matrix
      assert(rChebVal.rows() == this->mspSetup->bwdSize());
      assert(rChebVal.cols() == this->mspSetup->howmany());

      // Check for parity flipping operator
      int flip = this->flipsParity(integrator);

      // Do even/odd transforms
      for(int parity = 0; parity < 2; ++parity)
      {
         ParityTransformTools::extractParityModes(this->mTmpIn, physVal, this->parityBlocks(parity), physVal.rows());
         fftw_execute_r2r(this->fPlan((parity+flip)%2,parity), this->mTmpIn.data(), this->mTmpOut.data());
         ParityTransformTools::setParityModes(rChebVal, this->mTmpOut, this->parityBlocks(parity), physVal.rows());
      }

      if(integrator == SphereChebyshevFftwTransform::IntegratorType::INTG)
      {
         for(int parity = 0; parity < 2; ++parity)
         {
            // Rescale to remove FFT scaling
            ParityTransformTools::scaleParityModes(rChebVal, this->parityBlocks(parity), this->mspSetup->scale(), this->mspSetup->specSize());
         }

      } else
      {
         for(int parity = 0; parity < 2; ++parity)
         {
            ParityTransformTools::applyOperator(rChebVal, this->intgOp(integrator, parity), this->parityBlocks((parity+1)%2), this->mspSetup->scale(), this->mspSetup->specSize());
         }
      }

      #ifdef GEOMHDISCC_DEBUG
         rChebVal.bottomRows(this->mspSetup->padSize()).setConstant(std::numeric_limits<MHDFloat>::quiet_NaN());
      #endif //GEOMHDISCC_DEBUG
   }

   void SphereChebyshevFftwTransform::project(Matrix& rPhysVal, const Matrix& chebVal, SphereChebyshevFftwTransform::ProjectorType::Id projector, Arithmetics::Id arithId)
   {
      assert(arithId == Arithmetics::SET);

      // Assert that a mixed transform was not setup
      assert(this->mspSetup->type() == FftSetup::REAL);

      // assert on the padding size
      assert(this->mspSetup->padSize() >= 0);
      assert(this->mspSetup->bwdSize() - this->mspSetup->padSize() >= 0);

      // assert right sizes for input matrix
      assert(chebVal.rows() == this->mspSetup->bwdSize());
      assert(chebVal.cols() == this->mspSetup->howmany());

      // assert right sizes for output matrix
      assert(rPhysVal.rows() == this->mspSetup->fwdSize());
      assert(rPhysVal.cols() == this->mspSetup->howmany());

      // Check for parity flipping operator
      int flip = this->flipsParity(projector);

      // Loop over parity
      for(int parity = 0; parity < 2; ++parity)
      {
         // Compute first derivative
         if(projector == SphereChebyshevFftwTransform::ProjectorType::DIFF)
         {
            ParityTransformTools::extractParityModes(this->tmpInS(parity), chebVal, this->parityBlocks(parity), this->mspSetup->specSize());
            if((parity+1)%2 == 1)
            {
               this->tmpInS(parity).topRows(1).setZero();
            }
            this->tmpInS(parity).bottomRows(this->mspSetup->padSize()).setZero();
            Solver::internal::solveWrapper(this->tmpOutS(parity), this->solver(projector,(parity+1)%2), this->tmpInS(parity));
            this->mTmpIn.leftCols(this->tmpOutS(parity).cols()) = this->tmpOutS(parity);

         // Compute second derivative
         } else if(projector == SphereChebyshevFftwTransform::ProjectorType::DIFF2)
         {
            ParityTransformTools::extractParityModes(this->tmpInS(parity), chebVal, this->parityBlocks(parity), this->mspSetup->specSize());
            this->tmpInS(parity).topRows(1).setZero();
            this->tmpInS(parity).bottomRows(this->mspSetup->padSize()).setZero();
            Solver::internal::solveWrapper(this->tmpOutS(parity), this->solver(projector,parity), this->tmpInS(parity));
            this->mTmpIn.leftCols(this->tmpOutS(parity).cols()) = this->tmpOutS(parity);

         // Compute division by R
         } else if(projector == SphereChebyshevFftwTransform::ProjectorType::DIVR)
         {
            ParityTransformTools::extractParityModes(this->tmpInS(parity), chebVal, this->parityBlocks(parity), this->mspSetup->specSize());
            this->tmpInS(parity).bottomRows(this->mspSetup->padSize()).setZero();
            Solver::internal::solveWrapper(this->tmpOutS(parity), this->solver(projector,(parity+1)%2), this->tmpInS(parity));
            this->mTmpIn.leftCols(this->tmpOutS(parity).cols()) = this->tmpOutS(parity);

         // Compute division by R^2
         } else if(projector == SphereChebyshevFftwTransform::ProjectorType::DIVR2)
         {
            ParityTransformTools::extractParityModes(this->tmpInS(parity), chebVal, this->parityBlocks(parity), this->mspSetup->specSize());
            this->tmpInS(parity).bottomRows(this->mspSetup->padSize()).setZero();
            Solver::internal::solveWrapper(this->tmpOutS(parity), this->solver(projector,parity), this->tmpInS(parity));
            this->mTmpIn.leftCols(this->tmpOutS(parity).cols()) = this->tmpOutS(parity);

         // Compute D r projection
         } else if(projector == SphereChebyshevFftwTransform::ProjectorType::DIFFR)
         {
            ParityTransformTools::extractParityModes(this->tmpInS(parity), chebVal, this->parityBlocks(parity), this->mspSetup->specSize());
            this->tmpInS(parity) = this->solveOp(ProjectorType::DIVR, parity).leftCols(this->mspSetup->specSize())*this->tmpInS(parity).topRows(this->mspSetup->specSize());
            if(parity%2 == 1)
            {
               this->tmpInS(parity).topRows(1).setZero();
            }
            Solver::internal::solveWrapper(this->tmpOutS(parity), this->solver(ProjectorType::DIFF, parity), this->tmpInS(parity));
            this->mTmpIn.leftCols(this->tmpOutS(parity).cols()) = this->tmpOutS(parity);

         // Compute 1/r D r projection
         } else if(projector == SphereChebyshevFftwTransform::ProjectorType::DIVRDIFFR)
         {
            this->tmpInS(parity) = this->solveOp(ProjectorType::DIVR, parity).leftCols(this->mspSetup->specSize())*this->tmpInS(parity).topRows(this->mspSetup->specSize());
            if(parity == 1)
            {
               this->tmpInS(parity).topRows(1).setZero();
            }
            Solver::internal::solveWrapper(this->tmpOutS(parity), this->solver(ProjectorType::DIFF, parity), this->tmpInS(parity));
            Solver::internal::solveWrapper(this->tmpInS(parity), this->solver(ProjectorType::DIVR, (parity+1)%2), this->tmpOutS(parity));
            this->mTmpIn.leftCols(this->tmpInS(parity).cols()) = this->tmpInS(parity);

         // Compute radial laplacian projection
         } else if(projector == SphereChebyshevFftwTransform::ProjectorType::RADLAPL)
         {
            ParityTransformTools::extractParityModes(this->tmpInS(parity), chebVal, this->parityBlocks(parity), this->mspSetup->specSize());
            if((parity+1)%2 == 1)
            {
               this->tmpInS(parity).topRows(1).setZero();
            }
            this->tmpInS(parity).bottomRows(this->mspSetup->padSize()).setZero();
            Solver::internal::solveWrapper(this->tmpOutS(parity), this->solver(ProjectorType::DIFF,(parity+1)%2), this->tmpInS(parity));
            this->tmpInS(parity) = this->solveOp(ProjectorType::DIVR2, (parity+1)%2)*this->tmpOutS(parity);
            if(parity == 1)
            {
               this->tmpInS(parity).topRows(1).setZero();
            }
            Solver::internal::solveWrapper(this->tmpOutS(parity), this->solver(ProjectorType::DIFF, parity), this->tmpInS(parity));
            Solver::internal::solveWrapper(this->tmpInS(parity), this->solver(ProjectorType::DIVR2, parity), this->tmpOutS(parity));
            this->mTmpIn.leftCols(this->tmpInS(parity).cols()) = this->tmpInS(parity);

         // Compute simple projection
         } else
         {
            // Copy into other array
            ParityTransformTools::extractParityModes(this->mTmpIn, chebVal, this->parityBlocks(parity), this->mspSetup->specSize());
            this->mTmpIn.bottomRows(this->mspSetup->padSize()).setZero();
         }

         // Do transform
         fftw_execute_r2r(this->bPlan((parity+flip)%2, parity), this->mTmpIn.data(), this->mTmpOut.data());

         // Copy data into output variable
         ParityTransformTools::setParityModes(rPhysVal, this->mTmpOut, this->parityBlocks(parity), rPhysVal.rows());
      }
   }

   void SphereChebyshevFftwTransform::integrate(MatrixZ& rChebVal, const MatrixZ& physVal, SphereChebyshevFftwTransform::IntegratorType::Id integrator, Arithmetics::Id arithId)
   {
      assert(arithId == Arithmetics::SET);

      // Assert that a mixed transform was setup
      assert(this->mspSetup->type() == FftSetup::COMPONENT);

      // assert right sizes for input matrix
      assert(physVal.rows() == this->mspSetup->fwdSize());
      assert(physVal.cols() == this->mspSetup->howmany());

      // assert right sizes for output matrix
      assert(rChebVal.rows() == this->mspSetup->bwdSize());
      assert(rChebVal.cols() == this->mspSetup->howmany());

      // Check for parity flipping operator
      int flip = this->flipsParity(integrator);

      for(int component = 0; component < 2; ++component)
      {
         // Do even/odd transforms
         for(int parity = 0; parity < 2; ++parity)
         {
            ParityTransformTools::extractParityModes(this->mTmpIn, physVal, component, this->parityBlocks(parity), physVal.rows());
            fftw_execute_r2r(this->fPlan((parity+flip)%2,parity), this->mTmpIn.data(), this->mTmpOut.data());
            ParityTransformTools::setParityModes(rChebVal, this->mTmpOut, component, this->parityBlocks(parity), physVal.rows());
         }

         if(integrator == SphereChebyshevFftwTransform::IntegratorType::INTG)
         {
            for(int parity = 0; parity < 2; ++parity)
            {
               // Rescale to remove FFT scaling
               ParityTransformTools::scaleParityModes(rChebVal, component, this->parityBlocks(parity), this->mspSetup->scale(), physVal.rows());
               ParityTransformTools::extractParityModes(this->mTmpIn, rChebVal, component, this->parityBlocks(parity), physVal.rows());
               this->regularize(this->mTmpIn, parity);
               ParityTransformTools::setParityModes(rChebVal, this->mTmpIn, component, this->parityBlocks(parity), physVal.rows());
            }

         } else
         {
            for(int parity = 0; parity < 2; ++parity)
            {
               ParityTransformTools::applyOperator(rChebVal, component, this->intgOp(integrator,(parity+flip)%2), this->parityBlocks(parity), this->mspSetup->scale(), this->mspSetup->specSize());
            }
         }
      }

      #ifdef GEOMHDISCC_DEBUG
         rChebVal.bottomRows(this->mspSetup->padSize()).setConstant(std::numeric_limits<MHDFloat>::quiet_NaN());
      #endif //GEOMHDISCC_DEBUG
   }

   void SphereChebyshevFftwTransform::project(MatrixZ& rPhysVal, const MatrixZ& chebVal, SphereChebyshevFftwTransform::ProjectorType::Id projector, Arithmetics::Id arithId)
   {
      assert(arithId == Arithmetics::SET);

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

      // Check for parity flipping operator
      int flip = this->flipsParity(projector);

      for(int component = 0; component < 2; ++component)
      {
         // Loop over parity
         for(int parity = 0; parity < 2; ++parity)
         {
            // Compute first derivative
            if(projector == SphereChebyshevFftwTransform::ProjectorType::DIFF)
            {
               ParityTransformTools::extractParityModes(this->tmpInS(parity), chebVal, component, this->parityBlocks(parity), this->mspSetup->specSize());
               if((parity+1)%2 == 1)
               {
                  this->tmpInS(parity).topRows(1).setZero();
               }
               this->tmpInS(parity).bottomRows(this->mspSetup->padSize()).setZero();
               Solver::internal::solveWrapper(this->tmpOutS(parity), this->solver(projector,(parity+1)%2), this->tmpInS(parity));
               this->mTmpIn.leftCols(this->tmpOutS(parity).cols()) = this->tmpOutS(parity);

            // Compute second derivative
            } else if(projector == SphereChebyshevFftwTransform::ProjectorType::DIFF2)
            {
               ParityTransformTools::extractParityModes(this->tmpInS(parity), chebVal, component, this->parityBlocks(parity), this->mspSetup->specSize());
               this->tmpInS(parity).topRows(1).setZero();
               this->tmpInS(parity).bottomRows(this->mspSetup->padSize()).setZero();
               Solver::internal::solveWrapper(this->tmpOutS(parity), this->solver(projector,parity), this->tmpInS(parity));
               this->mTmpIn.leftCols(this->tmpOutS(parity).cols()) = this->tmpOutS(parity);

            // Compute division by R
            } else if(projector == SphereChebyshevFftwTransform::ProjectorType::DIVR)
            {
               ParityTransformTools::extractParityModes(this->tmpInS(parity), chebVal, component, this->parityBlocks(parity), this->mspSetup->specSize());
               this->tmpInS(parity).bottomRows(this->mspSetup->padSize()).setZero();

               // With exact r^l value should be zero, compute and correct regularity error (without it it will generate a flat spectrum for c/r)
               if(parity == 0)
               {
                  ParityTransformTools::correctRegularity(this->tmpInS(parity), this->mspSetup->specSize());
               }
               Solver::internal::solveWrapper(this->tmpOutS(parity), this->solver(projector,(parity+1)%2), this->tmpInS(parity));
               this->mTmpIn.leftCols(this->tmpOutS(parity).cols()) = this->tmpOutS(parity);

               // Compute division by R^2
            } else if(projector == SphereChebyshevFftwTransform::ProjectorType::DIVR2)
            {
               ParityTransformTools::extractParityModes(this->tmpInS(parity), chebVal, component, this->parityBlocks(parity), this->mspSetup->specSize());
               this->tmpInS(parity).bottomRows(this->mspSetup->padSize()).setZero();
               Solver::internal::solveWrapper(this->tmpOutS(parity), this->solver(projector,parity), this->tmpInS(parity));
               this->mTmpIn.leftCols(this->tmpOutS(parity).cols()) = this->tmpOutS(parity);

               // Compute D r projection
            } else if(projector == SphereChebyshevFftwTransform::ProjectorType::DIFFR)
            {
               ParityTransformTools::extractParityModes(this->tmpInS(parity), chebVal, component, this->parityBlocks(parity), this->mspSetup->specSize());
               this->tmpInS(parity) = this->solveOp(ProjectorType::DIVR, parity).leftCols(this->mspSetup->specSize())*this->tmpInS(parity).topRows(this->mspSetup->specSize());
               if(parity%2 == 1)
               {
                  this->tmpInS(parity).topRows(1).setZero();
               }
               Solver::internal::solveWrapper(this->tmpOutS(parity), this->solver(ProjectorType::DIFF, parity), this->tmpInS(parity));
               this->mTmpIn.leftCols(this->tmpOutS(parity).cols()) = this->tmpOutS(parity);

               // Compute 1/r D r projection
            } else if(projector == SphereChebyshevFftwTransform::ProjectorType::DIVRDIFFR)
            {
               ParityTransformTools::extractParityModes(this->tmpInS(parity), chebVal, component, this->parityBlocks(parity), this->mspSetup->specSize());
               this->tmpInS(parity) = this->solveOp(ProjectorType::DIVR, parity).leftCols(this->mspSetup->specSize())*this->tmpInS(parity).topRows(this->mspSetup->specSize());
               if(parity == 1)
               {
                  this->tmpInS(parity).topRows(1).setZero();
               }
               Solver::internal::solveWrapper(this->tmpOutS(parity), this->solver(ProjectorType::DIFF, parity), this->tmpInS(parity));
               Solver::internal::solveWrapper(this->tmpInS(parity), this->solver(ProjectorType::DIVR, (parity+1)%2), this->tmpOutS(parity));
               this->mTmpIn.leftCols(this->tmpInS(parity).cols()) = this->tmpInS(parity);

               // Compute radial laplacian projection
            } else if(projector == SphereChebyshevFftwTransform::ProjectorType::RADLAPL)
            {
               ParityTransformTools::extractParityModes(this->tmpInS(parity), chebVal, component, this->parityBlocks(parity), this->mspSetup->specSize());
               if((parity+1)%2 == 1)
               {
                  this->tmpInS(parity).topRows(1).setZero();
               }
               this->tmpInS(parity).bottomRows(this->mspSetup->padSize()).setZero();
               Solver::internal::solveWrapper(this->tmpOutS(parity), this->solver(ProjectorType::DIFF,(parity+1)%2), this->tmpInS(parity));
               this->tmpInS(parity) = this->solveOp(ProjectorType::DIVR2, (parity+1)%2)*this->tmpOutS(parity);
               if(parity == 1)
               {
                  this->tmpInS(parity).topRows(1).setZero();
               }
               Solver::internal::solveWrapper(this->tmpOutS(parity), this->solver(ProjectorType::DIFF, parity), this->tmpInS(parity));
               Solver::internal::solveWrapper(this->tmpInS(parity), this->solver(ProjectorType::DIVR2, parity), this->tmpOutS(parity));
               this->mTmpIn.leftCols(this->tmpInS(parity).cols()) = this->tmpInS(parity);

            // Compute simple projection
            } else
            {
               // Copy into other array
               ParityTransformTools::extractParityModes(this->mTmpIn, chebVal, component, this->parityBlocks(parity), this->mspSetup->specSize());
               this->mTmpIn.bottomRows(this->mspSetup->padSize()).setZero();
            }

            // Do transform
            fftw_execute_r2r(this->bPlan((parity+flip)%2, parity), this->mTmpIn.data(), this->mTmpOut.data());

            // Copy data into output variable
            ParityTransformTools::setParityModes(rPhysVal, this->mTmpOut, component, this->parityBlocks(parity), rPhysVal.rows());
         }
      }
   }

   void SphereChebyshevFftwTransform::integrate_full(Matrix& rChebVal, const Matrix& physVal, SphereChebyshevFftwTransform::IntegratorType::Id integrator, Arithmetics::Id arithId)
   {
      assert(arithId == Arithmetics::SET);

      // Assert that a mixed transform was not setup
      assert(this->mspSetup->type() == FftSetup::REAL);

      // assert right sizes for input matrix
      assert(physVal.rows() == this->mspSetup->fwdSize());
      assert(physVal.cols() == this->mspSetup->howmany());

      // assert right sizes for output matrix
      assert(rChebVal.rows() == this->mspSetup->bwdSize());
      assert(rChebVal.cols() == this->mspSetup->howmany());

      // Check for parity flipping operator
      int flip = this->flipsParity(integrator);

      // Do even/odd transforms
      for(int parity = 0; parity < 2; ++parity)
      {
         ParityTransformTools::extractParityModes(this->mTmpIn, physVal, this->parityBlocks(parity), physVal.rows());
         fftw_execute_r2r(this->fPlan((parity+flip)%2,parity), this->mTmpIn.data(), this->mTmpOut.data());
         ParityTransformTools::setParityModes(rChebVal, this->mTmpOut, this->parityBlocks(parity), physVal.rows());
      }

      if(integrator == SphereChebyshevFftwTransform::IntegratorType::INTG)
      {
         for(int parity = 0; parity < 2; ++parity)
         {
            // Rescale to remove FFT scaling
            ParityTransformTools::scaleParityModes(rChebVal, this->parityBlocks(parity), this->mspSetup->scale(), rChebVal.rows());
         }

      } else
      {
         for(int parity = 0; parity < 2; ++parity)
         {
            ParityTransformTools::applyOperator(rChebVal, this->intgOp(integrator, parity), this->parityBlocks(parity), this->mspSetup->scale(), rChebVal.rows());
         }
      }
   }

   void SphereChebyshevFftwTransform::integrate_full(MatrixZ& rChebVal, const MatrixZ& physVal, SphereChebyshevFftwTransform::IntegratorType::Id integrator, Arithmetics::Id arithId)
   {
      assert(arithId == Arithmetics::SET);

      // Assert that a mixed transform was setup
      assert(this->mspSetup->type() == FftSetup::COMPONENT);

      // assert right sizes for input matrix
      assert(physVal.rows() == this->mspSetup->fwdSize());
      assert(physVal.cols() == this->mspSetup->howmany());

      // assert right sizes for output matrix
      assert(rChebVal.rows() == this->mspSetup->bwdSize());
      assert(rChebVal.cols() == this->mspSetup->howmany());

      // Check for parity flipping operator
      int flip = this->flipsParity(integrator);

      for(int component = 0; component < 2; ++component)
      {
         // Do even/odd transforms
         for(int parity = 0; parity < 2; ++parity)
         {
            ParityTransformTools::extractParityModes(this->mTmpIn, physVal, component, this->parityBlocks(parity), physVal.rows());
            fftw_execute_r2r(this->fPlan((parity+flip)%2,parity), this->mTmpIn.data(), this->mTmpOut.data());
            ParityTransformTools::setParityModes(rChebVal, this->mTmpOut, component, this->parityBlocks(parity), physVal.rows());
         }

         if(integrator == SphereChebyshevFftwTransform::IntegratorType::INTG)
         {
            for(int parity = 0; parity < 2; ++parity)
            {
               // Rescale to remove FFT scaling
               ParityTransformTools::scaleParityModes(rChebVal, component, this->parityBlocks(parity), this->mspSetup->scale(), rChebVal.rows());
            }

         } else
         {
            for(int parity = 0; parity < 2; ++parity)
            {
               ParityTransformTools::applyOperator(rChebVal, component, this->intgOp(integrator,parity), this->parityBlocks(parity), this->mspSetup->scale(), rChebVal.rows());
            }
         }
      }
   }

   void SphereChebyshevFftwTransform::integrate_energy(Matrix& rChebVal, const Matrix& physVal, SphereChebyshevFftwTransform::IntegratorType::Id integrator, Arithmetics::Id arithId)
   {
      assert(arithId == Arithmetics::SET);

      // Assert that a mixed transform was not setup
      assert(this->mspSetup->type() == FftSetup::REAL);

      // assert right sizes for input matrix
      assert(physVal.rows() == this->mspSetup->fwdSize());
      assert(physVal.cols() == this->mspSetup->howmany());

      // assert right sizes for output matrix
      assert(rChebVal.rows() == this->mspSetup->bwdSize());
      assert(rChebVal.cols() == this->mspSetup->howmany());

      // Do even/odd transforms
      for(int parity = 0; parity < 2; ++parity)
      {
         ParityTransformTools::extractParityModes(this->mTmpIn, physVal, this->parityBlocks(parity), physVal.rows());
         fftw_execute_r2r(this->fPlan(0,0), this->mTmpIn.data(), this->mTmpOut.data());
         ParityTransformTools::setParityModes(rChebVal, this->mTmpOut, this->parityBlocks(parity), physVal.rows());
         // Rescale to remove FFT scaling
         ParityTransformTools::scaleParityModes(rChebVal, this->parityBlocks(parity), this->mspSetup->scale(), rChebVal.rows());
      }
   }

   void SphereChebyshevFftwTransform::integrate_energy(MatrixZ& rChebVal, const MatrixZ& physVal, SphereChebyshevFftwTransform::IntegratorType::Id integrator, Arithmetics::Id arithId)
   {
      assert(arithId == Arithmetics::SET);

      // Assert that a mixed transform was setup
      assert(this->mspSetup->type() == FftSetup::COMPONENT);

      // assert right sizes for input matrix
      assert(physVal.rows() == this->mspSetup->fwdSize());
      assert(physVal.cols() == this->mspSetup->howmany());

      // assert right sizes for output matrix
      assert(rChebVal.rows() == this->mspSetup->bwdSize());
      assert(rChebVal.cols() == this->mspSetup->howmany());

      // Energy transform is only even
      for(int parity = 0; parity < 2; ++parity)
      {
         ParityTransformTools::extractParityModes(this->mTmpIn, physVal, 1, this->parityBlocks(parity), physVal.rows());
         fftw_execute_r2r(this->fPlan(0,0), this->mTmpIn.data(), this->mTmpOut.data());
         ParityTransformTools::setParityModes(rChebVal, this->mTmpOut, 1, this->parityBlocks(parity), physVal.rows());
         // Rescale to remove FFT scaling
         ParityTransformTools::scaleParityModes(rChebVal, 1, this->parityBlocks(parity), this->mspSetup->scale(), rChebVal.rows());
      }
      rChebVal.imag().setZero();
   }

   void SphereChebyshevFftwTransform::regularize(Matrix& rData, const int parity)
   {
      SphereChebyshevFftwTransform::RegularityType::Id reg = SphereChebyshevFftwTransform::RegularityType::REG2;

      Matrix mRegIn = rData.bottomRows(rData.rows() - this->regularitySize(reg, parity));
      Matrix mRegOut = Matrix::Zero(rData.rows() - this->regularitySize(reg, parity), rData.cols());

      Solver::internal::solveWrapper(mRegOut, this->regSolver(reg, parity), mRegIn);

      rData.topRows(this->mspSetup->specSize()) = this->regOp(reg, parity)*mRegOut.topRows(this->mspSetup->specSize() - this->regularitySize(reg, parity));
   }

#ifdef GEOMHDISCC_STORAGEPROFILE
   MHDFloat SphereChebyshevFftwTransform::requiredStorage() const
   {
      MHDFloat mem = 0.0;

      // Storage required for the fftw plans 
      mem += 8.0*2.0;

      return mem;
   }
#endif // GEOMHDISCC_STORAGEPROFILE

}
}
