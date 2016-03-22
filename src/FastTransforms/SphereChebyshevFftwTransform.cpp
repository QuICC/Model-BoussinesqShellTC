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
      // Abort if the imposed regularity is too small
      if(RegularityType::REGMAX < RegularityType::REG3)
      {
         throw Exception("Imposed regularity is too small for toroidal/poloidal operators!");
      }

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
      this->mRegularityShift.insert(std::make_pair(IntegratorType::INTG, std::make_pair(0,0)));
      // Multiplication by R
      this->mIntgOp.insert(std::make_pair(IntegratorType::INTGR, opPair));
      this->mIntegratorFlips.insert(std::make_pair(IntegratorType::INTGR, 1));
      this->mRegularityShift.insert(std::make_pair(IntegratorType::INTGR, std::make_pair(0,0)));
      // QST Q operator (4th order)
      this->mIntgOp.insert(std::make_pair(IntegratorType::INTGQ4, opPair));
      this->mIntegratorFlips.insert(std::make_pair(IntegratorType::INTGQ4, 1));
      this->mRegularityShift.insert(std::make_pair(IntegratorType::INTGQ4, std::make_pair(4,4)));
      // QST S operator (4th order)
      this->mIntgOp.insert(std::make_pair(IntegratorType::INTGS4, opPair));
      this->mIntegratorFlips.insert(std::make_pair(IntegratorType::INTGS4, 1));
      this->mRegularityShift.insert(std::make_pair(IntegratorType::INTGS4, std::make_pair(4,4)));
      // QST T operator
      this->mIntgOp.insert(std::make_pair(IntegratorType::INTGT, opPair));
      this->mIntegratorFlips.insert(std::make_pair(IntegratorType::INTGT, 0));
      this->mRegularityShift.insert(std::make_pair(IntegratorType::INTGT, std::make_pair(2,2)));
      // QST Q operator (2th order)
      this->mIntgOp.insert(std::make_pair(IntegratorType::INTGQ2, opPair));
      this->mIntegratorFlips.insert(std::make_pair(IntegratorType::INTGQ2, 1));
      this->mRegularityShift.insert(std::make_pair(IntegratorType::INTGQ2, std::make_pair(2,2)));
      // QST S operator (2th order)
      this->mIntgOp.insert(std::make_pair(IntegratorType::INTGS2, opPair));
      this->mIntegratorFlips.insert(std::make_pair(IntegratorType::INTGS2, 1));
      this->mRegularityShift.insert(std::make_pair(IntegratorType::INTGS2, std::make_pair(2,2)));

      //
      // Initialise regularity operators
      //
      for(int i = 0; i <= static_cast<int>(RegularityType::REGMAX); ++i)
      {
         this->mRegOp.insert(std::make_pair(static_cast<RegularityType::Id>(i), opPair));
         this->mRegularitySize.insert(std::make_pair(static_cast<RegularityType::Id>(i), std::make_pair(1 + i/2, i/2 + i%2)));
      }

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
      for(int i = 0; i <= static_cast<int>(RegularityType::REGMAX); ++i)
      {
         this->mRegSolveOp.insert(std::make_pair(static_cast<RegularityType::Id>(i), opPair));
      }

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
      for(int i = 0; i <= static_cast<int>(RegularityType::REGMAX); ++i)
      {
         pSolE = SharedPtrMacro<Solver::SparseSelector<SparseMatrix>::Type>(new Solver::SparseSelector<SparseMatrix>::Type());
         pSolO = SharedPtrMacro<Solver::SparseSelector<SparseMatrix>::Type>(new Solver::SparseSelector<SparseMatrix>::Type());
         this->mRegSolver.insert(std::make_pair(static_cast<RegularityType::Id>(i), std::make_pair(pSolE,pSolO)));
      }

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

         // Prepare regularity operators
         pValue = PyLong_FromLong(this->mspSetup->specSize());
         PyTuple_SetItem(pRegArgs, 0, pValue);
         PyTuple_SetItem(pRegArgs, 1, PyLong_FromLong(parity));
         pValue = PyDict_New();
         PyDict_SetItem(pValue, PyLong_FromLong(0), PyLong_FromLong(0));
         PyTuple_SetItem(pRegArgs, 2, pValue);
         PyTuple_SetItem(pRegArgs, 3, PyLong_FromLong(0));

         // Call stencil for regularity solver
         std::vector<int> regBc;
         regBc.push_back(-14); // Boundary condition for 0th order regularity
         regBc.push_back(-15); // Boundary condition for 1st order regularity
         regBc.push_back(-22); // Boundary condition for 2nd order regularity
         regBc.push_back(-23); // Boundary condition for 3rd order regularity
         regBc.push_back(-30); // Boundary condition for 4th order regularity
         regBc.push_back(-31); // Boundary condition for 5th order regularity
         regBc.push_back(-40); // Boundary condition for 6th order regularity
         regBc.push_back(-41); // Boundary condition for 7th order regularity
         regBc.push_back(-50); // Boundary condition for 8th order regularity
         for(int i = 0; i <= static_cast<int>(RegularityType::REGMAX); ++i)
         {
            if(this->regularitySize(static_cast<RegularityType::Id>(i), parity) > 0)
            {
               pValue = PyTuple_GetItem(pRegArgs, 2);
               PyDict_SetItem(pValue, PyLong_FromLong(0), PyLong_FromLong(regBc.at(i)));
               PythonWrapper::setFunction("stencil");
               pValue = PythonWrapper::callFunction(pRegArgs);
               // Fill matrix
               PythonWrapper::fillMatrix(this->regOp(static_cast<RegularityType::Id>(i), parity), pValue);
               Py_DECREF(pValue);
            }
         }

         // Prepare regularity solvers
         pValue = PyLong_FromLong(this->mspSetup->fwdSize());
         PyTuple_SetItem(pRegArgs, 0, pValue);
         PyTuple_SetItem(pRegArgs, 3, PyLong_FromLong(1));

         // Call stencil for regularity solver
         for(int i = 0; i <= static_cast<int>(RegularityType::REGMAX); ++i)
         {
            if(this->regularitySize(static_cast<RegularityType::Id>(i), parity) > 0)
            {
               pValue = PyTuple_GetItem(pRegArgs, 2);
               PyDict_SetItem(pValue, PyLong_FromLong(0), PyLong_FromLong(regBc.at(i)));
               PyDict_SetItem(pValue, PyUnicode_FromString("rt"), PyLong_FromLong(this->regularitySize(static_cast<RegularityType::Id>(i), parity)));
               PythonWrapper::setFunction("stencil");
               pValue = PythonWrapper::callFunction(pRegArgs);
               // Fill matrix
               PythonWrapper::fillMatrix(this->regSolveOp(static_cast<RegularityType::Id>(i), parity), pValue);
               Py_DECREF(pValue);
            }
         }

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

         // Initialize solver and factorize regularity stencil
         for(int i = 0; i <= static_cast<int>(RegularityType::REGMAX); ++i)
         {
            if(this->regularitySize(static_cast<RegularityType::Id>(i), parity) > 0)
            {
               this->regSolver(static_cast<RegularityType::Id>(i), parity).compute(this->regSolveOp(static_cast<RegularityType::Id>(i), parity));
               // Check for successful factorisation
               if(this->regSolver(static_cast<RegularityType::Id>(i), parity).info() != Eigen::Success)
               {
                  throw Exception("Factorization of regularity stencil failed!");
               }
            }
         }
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
            ParityTransformTools::scaleParityModes(rChebVal, this->parityBlocks(parity), this->mspSetup->scale(), rChebVal.rows());
         }

      } else
      {
         for(int parity = 0; parity < 2; ++parity)
         {
            ParityTransformTools::applyOperator(rChebVal, this->intgOp(integrator, parity), this->parityBlocks((parity+1)%2), this->mspSetup->scale(), rChebVal.rows());
         }
      }

      // Regularize solution
      for(int parity = 0; parity < 2; ++parity)
      {
         ParityTransformTools::extractParityModes(this->mTmpIn, rChebVal, this->parityBlocks(parity), rChebVal.rows());
         this->regularize(this->mTmpIn, parity, this->regularityShift(integrator, parity));
         ParityTransformTools::setParityModes(rChebVal, this->mTmpIn, this->parityBlocks(parity), rChebVal.rows());
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
         ParityTransformTools::extractParityModes(this->tmpInS(parity), chebVal, this->parityBlocks(parity), this->mspSetup->specSize());
         this->tmpInS(parity).bottomRows(this->mspSetup->padSize()).setZero();
         this->regularize(this->tmpInS(parity), parity, 0);

         // Compute first derivative
         if(projector == SphereChebyshevFftwTransform::ProjectorType::DIFF)
         {
            if((parity+1)%2 == 1)
            {
               this->tmpInS(parity).topRows(1).setZero();
            }
            Solver::internal::solveWrapper(this->tmpOutS(parity), this->solver(projector,(parity+1)%2), this->tmpInS(parity));
            this->mTmpIn.leftCols(this->tmpOutS(parity).cols()) = this->tmpOutS(parity);

         // Compute second derivative
         } else if(projector == SphereChebyshevFftwTransform::ProjectorType::DIFF2)
         {
            this->tmpInS(parity).topRows(1).setZero();
            Solver::internal::solveWrapper(this->tmpOutS(parity), this->solver(projector,parity), this->tmpInS(parity));
            this->mTmpIn.leftCols(this->tmpOutS(parity).cols()) = this->tmpOutS(parity);

         // Compute division by R
         } else if(projector == SphereChebyshevFftwTransform::ProjectorType::DIVR)
         {
            Solver::internal::solveWrapper(this->tmpOutS(parity), this->solver(projector,(parity+1)%2), this->tmpInS(parity));
            this->mTmpIn.leftCols(this->tmpOutS(parity).cols()) = this->tmpOutS(parity);

         // Compute division by R^2
         } else if(projector == SphereChebyshevFftwTransform::ProjectorType::DIVR2)
         {
            Solver::internal::solveWrapper(this->tmpOutS(parity), this->solver(projector,parity), this->tmpInS(parity));
            this->mTmpIn.leftCols(this->tmpOutS(parity).cols()) = this->tmpOutS(parity);

         // Compute D r projection
         } else if(projector == SphereChebyshevFftwTransform::ProjectorType::DIFFR)
         {
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
            if((parity+1)%2 == 1)
            {
               this->tmpInS(parity).topRows(1).setZero();
            }
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
            this->mTmpIn.leftCols(this->tmpInS(parity).cols()) = this->tmpInS(parity);
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
               ParityTransformTools::scaleParityModes(rChebVal, component, this->parityBlocks(parity), this->mspSetup->scale(), rChebVal.rows());
            }

         } else
         {
            for(int parity = 0; parity < 2; ++parity)
            {
               ParityTransformTools::applyOperator(rChebVal, component, this->intgOp(integrator,(parity+flip)%2), this->parityBlocks(parity), this->mspSetup->scale(), rChebVal.rows());
            }
         }

         // Regularize solution
         for(int parity = 0; parity < 2; ++parity)
         {
            ParityTransformTools::extractParityModes(this->mTmpIn, rChebVal, component, this->parityBlocks(parity), rChebVal.rows());
            this->regularize(this->mTmpIn, parity, this->regularityShift(integrator, parity));
            ParityTransformTools::setParityModes(rChebVal, this->mTmpIn, component, this->parityBlocks(parity), rChebVal.rows());
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
            ParityTransformTools::extractParityModes(this->tmpInS(parity), chebVal, component, this->parityBlocks(parity), this->mspSetup->specSize());
            this->tmpInS(parity).bottomRows(this->mspSetup->padSize()).setZero();
            this->regularize(this->tmpInS(parity), parity, 0);

            // Compute first derivative
            if(projector == SphereChebyshevFftwTransform::ProjectorType::DIFF)
            {
               if((parity+1)%2 == 1)
               {
                  this->tmpInS(parity).topRows(1).setZero();
               }
               Solver::internal::solveWrapper(this->tmpOutS(parity), this->solver(projector,(parity+1)%2), this->tmpInS(parity));
               this->mTmpIn.leftCols(this->tmpOutS(parity).cols()) = this->tmpOutS(parity);

            // Compute second derivative
            } else if(projector == SphereChebyshevFftwTransform::ProjectorType::DIFF2)
            {
               this->tmpInS(parity).topRows(1).setZero();
               Solver::internal::solveWrapper(this->tmpOutS(parity), this->solver(projector,parity), this->tmpInS(parity));
               this->mTmpIn.leftCols(this->tmpOutS(parity).cols()) = this->tmpOutS(parity);

            // Compute division by R
            } else if(projector == SphereChebyshevFftwTransform::ProjectorType::DIVR)
            {
               Solver::internal::solveWrapper(this->tmpOutS(parity), this->solver(projector,(parity+1)%2), this->tmpInS(parity));
               this->mTmpIn.leftCols(this->tmpOutS(parity).cols()) = this->tmpOutS(parity);

               // Compute division by R^2
            } else if(projector == SphereChebyshevFftwTransform::ProjectorType::DIVR2)
            {
               Solver::internal::solveWrapper(this->tmpOutS(parity), this->solver(projector,parity), this->tmpInS(parity));
               this->mTmpIn.leftCols(this->tmpOutS(parity).cols()) = this->tmpOutS(parity);

               // Compute D r projection
            } else if(projector == SphereChebyshevFftwTransform::ProjectorType::DIFFR)
            {
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
               if((parity+1)%2 == 1)
               {
                  this->tmpInS(parity).topRows(1).setZero();
               }
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
               this->mTmpIn.leftCols(this->tmpInS(parity).cols()) = this->tmpInS(parity);
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

   void SphereChebyshevFftwTransform::regularizeBlock(Matrix& rData, const int start, const int cols, const int parity, const SphereChebyshevFftwTransform::RegularityType::Id reg)
   {
      Matrix regIn = rData.block(this->regularitySize(reg, parity), start, rData.rows() - this->regularitySize(reg, parity), cols);
      Matrix regOut = Matrix::Zero(rData.rows() - this->regularitySize(reg, parity), cols);

      Solver::internal::solveWrapper(regOut, this->regSolver(reg, parity), regIn);

      rData.block(0, start, this->mspSetup->specSize(), cols) = this->regOp(reg, parity)*regOut.topRows(this->mspSetup->specSize() - this->regularitySize(reg, parity));
   }

   void SphereChebyshevFftwTransform::regularize(Matrix& rData, const int parity, const int regShift)
   {
      int start = 0;
      int cols = rData.cols();

      // Even special cases
      int idx = 0;
      // Special case for l = 0, or l = 1
      if(this->parityBlocks(parity)(idx, 2) == parity && regShift == 0)
      {
         start += this->parityBlocks(parity)(idx,1);
         cols -= this->parityBlocks(parity)(idx,1);
         idx++;
      }

      bool hasSpecial = (cols > 0);
      while(hasSpecial)
      {
         // Solution is less regular than max regularity
         if(this->parityBlocks(parity)(idx,2) - 1 + regShift < static_cast<int>(SphereChebyshevFftwTransform::RegularityType::REGMAX))
         {
            this->regularizeBlock(rData, start, this->parityBlocks(parity)(idx,1), parity, static_cast<SphereChebyshevFftwTransform::RegularityType::Id>(this->parityBlocks(parity)(idx,2) - 1 + regShift));
            start += this->parityBlocks(parity)(idx,1);
            cols -= this->parityBlocks(parity)(idx,1);
            idx++;

            if(cols < 1)
            {
               hasSpecial = false;
            }

         // Special cases are done
         } else
         {
            hasSpecial = false;
         }
      }

      if(cols > 0)
      {
         this->regularizeBlock(rData, start, cols, parity, SphereChebyshevFftwTransform::RegularityType::REGMAX);
      }
   }

   void SphereChebyshevFftwTransform::checkRegularity(const Matrix& data, const int rows)
   {
//      Array err = Array::Zero(data.cols());
//      for(int i = rows-1; i > 0; --i)
//      {
//         err.transpose() += std::pow(-1,i+1)*data.row(i);
//      }
//      err *= 2.0;
//      std::cerr << "l = 0: " << err(0) - data(0,0) << " l> 0: " << (err.bottomRows(err.rows()-1).transpose() - data.row(0).rightCols(data.cols()-1)).array().abs().maxCoeff() << std::endl;
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
