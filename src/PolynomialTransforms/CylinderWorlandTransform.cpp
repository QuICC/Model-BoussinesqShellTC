/** 
 * @file CylinderWorlandTransform.cpp
 * @brief Source of the implementation of the Worland transform in a cylinder
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "PolynomialTransforms/CylinderWorlandTransform.hpp"

// Project includes
//
#include "Exceptions/Exception.hpp"
#include "Base/MathConstants.hpp"
#include "Quadratures/WorlandChebyshevRule.hpp"
#include "PolynomialTransforms/WorlandPolynomial.hpp"
#include "Python/PythonWrapper.hpp"

namespace GeoMHDiSCC {

namespace Transform {

   Array CylinderWorlandTransform::generateGrid(const int size)
   {
      // Initialise grid storage
      Array grid(size);
      Array tmp(size);

      WorlandChebyshevRule::computeQuadrature(grid, tmp, size);

      return grid;
   }

   CylinderWorlandTransform::CylinderWorlandTransform()
   {
      // Initialise the Python interpreter wrapper
      PythonWrapper::init();
   }

   CylinderWorlandTransform::~CylinderWorlandTransform()
   {
   }

   void CylinderWorlandTransform::init(CylinderWorlandTransform::SharedSetupType spSetup)
   {
      // Store the shared pointer to setup object
      this->mspSetup = spSetup;

      // Initialise the quadrature grid and weights and operators
      this->initOperators();
   }

   void CylinderWorlandTransform::requiredOptions(std::set<NonDimensional::Id>& list, const Dimensions::Transform::Id dimId) const
   {
      //
      // No possible options
      //
   }

   void CylinderWorlandTransform::setOptions(const std::map<NonDimensional::Id, MHDFloat>& options, const Dimensions::Transform::Id dimId)
   {
      //
      // No possible options
      //
   }

   const Array& CylinderWorlandTransform::meshGrid() const
   {
      if(this->mGrid.size() == 0 || this->mWeights.size() == 0)
      {
         throw Exception("CylinderWorlandTransform has not been initialised!");
      }

      return this->mGrid;
   }

   void CylinderWorlandTransform::initOperators()
   {
      // Initialise python wrapper
      PythonWrapper::import("quicc.geometry.cylindrical.cylinder_radius_worland");

      this->mGrid.resize(this->mspSetup->fwdSize());
      this->mWeights.resize(this->mspSetup->fwdSize());

      // Set the grid and weights
      internal::Array igrid, iweights;
      WorlandChebyshevRule::computeQuadrature(this->mGrid, this->mWeights, igrid, iweights, this->mspSetup->fwdSize());

      // Reserve storage for the projectors
      this->mProjOp.insert(std::make_pair(ProjectorType::PROJ,std::vector<Matrix>()));
      this->mProjOp.find(ProjectorType::PROJ)->second.reserve(this->mspSetup->slow().size());
      this->mProjOp.insert(std::make_pair(ProjectorType::DIVRM0,std::vector<Matrix>()));
      this->mProjOp.find(ProjectorType::DIVRM0)->second.reserve(this->mspSetup->slow().size());
      this->mProjOp.insert(std::make_pair(ProjectorType::DIFF,std::vector<Matrix>()));
      this->mProjOp.find(ProjectorType::DIFF)->second.reserve(this->mspSetup->slow().size());
      this->mProjOp.insert(std::make_pair(ProjectorType::DIFF_PROJ,std::vector<Matrix>()));
      this->mProjOp.find(ProjectorType::DIFF_PROJ)->second.reserve(this->mspSetup->slow().size());
      this->mProjOp.insert(std::make_pair(ProjectorType::DIVRDIFFR,std::vector<Matrix>()));
      this->mProjOp.find(ProjectorType::DIVRDIFFR)->second.reserve(this->mspSetup->slow().size());
      this->mProjOp.insert(std::make_pair(ProjectorType::LAPLH,std::vector<Matrix>()));
      this->mProjOp.find(ProjectorType::LAPLH)->second.reserve(this->mspSetup->slow().size());
      this->mProjOp.insert(std::make_pair(ProjectorType::LAPLH_DIVRDIFFR,std::vector<Matrix>()));
      this->mProjOp.find(ProjectorType::LAPLH_DIVRDIFFR)->second.reserve(this->mspSetup->slow().size());
      this->mProjOp.insert(std::make_pair(ProjectorType::DIVRLAPLHM0,std::vector<Matrix>()));
      this->mProjOp.find(ProjectorType::DIVRLAPLHM0)->second.reserve(this->mspSetup->slow().size());
      this->mProjOp.insert(std::make_pair(ProjectorType::DIFFLAPLH,std::vector<Matrix>()));
      this->mProjOp.find(ProjectorType::DIFFLAPLH)->second.reserve(this->mspSetup->slow().size());
      this->mProjOp.insert(std::make_pair(ProjectorType::DIFFLAPLH_DIFFDIVRDIFFR,std::vector<Matrix>()));
      this->mProjOp.find(ProjectorType::DIFFLAPLH_DIFFDIVRDIFFR)->second.reserve(this->mspSetup->slow().size());

      // Reserve storage for the weighted projectors 
      this->mIntgOp.insert(std::make_pair(IntegratorType::INTG,std::vector<Matrix>()));
      this->mIntgOp.find(IntegratorType::INTG)->second.reserve(this->mspSetup->slow().size());
      this->mIntgOp.insert(std::make_pair(IntegratorType::INTGI4DIVRM0,std::vector<Matrix>()));
      this->mIntgOp.find(IntegratorType::INTGI4DIVRM0)->second.reserve(this->mspSetup->slow().size());
      this->mIntgOp.insert(std::make_pair(IntegratorType::INTGI4DIVRDIFFR_I2,std::vector<Matrix>()));
      this->mIntgOp.find(IntegratorType::INTGI4DIVRDIFFR_I2)->second.reserve(this->mspSetup->slow().size());
      this->mIntgOp.insert(std::make_pair(IntegratorType::INTGI6DIVRM0,std::vector<Matrix>()));
      this->mIntgOp.find(IntegratorType::INTGI6DIVRM0)->second.reserve(this->mspSetup->slow().size());
      this->mIntgOp.insert(std::make_pair(IntegratorType::INTGI6DIVRDIFFR_I4,std::vector<Matrix>()));
      this->mIntgOp.find(IntegratorType::INTGI6DIVRDIFFR_I4)->second.reserve(this->mspSetup->slow().size());
      this->mIntgOp.insert(std::make_pair(IntegratorType::INTGI6LAPLH_I4DR,std::vector<Matrix>()));
      this->mIntgOp.find(IntegratorType::INTGI6LAPLH_I4DR)->second.reserve(this->mspSetup->slow().size());
      this->mIntgOp.insert(std::make_pair(IntegratorType::INTGINVLAPLH,std::vector<Matrix>()));
      this->mIntgOp.find(IntegratorType::INTGINVLAPLH)->second.reserve(this->mspSetup->slow().size());

      // Reserve storage for the solver operators 
      this->mSolveOp.insert(std::make_pair(ProjectorType::INVLAPLH,std::vector<SparseMatrix>()));
      this->mSolveOp.find(ProjectorType::INVLAPLH)->second.reserve(this->mspSetup->slow().size());
      this->mTriSolver.insert(std::make_pair(ProjectorType::INVLAPLH,std::vector<SharedPtrMacro<Solver::SparseTriSelector<SparseMatrix>::Type> >()));
      this->mTriSolver.find(ProjectorType::INVLAPLH)->second.reserve(this->mspSetup->slow().size());

      // Prepare arguments to Python matrices call
      PyObject *pArgs, *pValue;
      pArgs = PyTuple_New(4);

      // Loop over harmonic degrees
      Matrix op;
      internal::Matrix  itmp, ipoly;
      for(int iL = 0; iL < this->mspSetup->slow().size(); iL++)
      {
         int l = this->mspSetup->slow()(iL);

         // Allocate memory for the projectors
         this->mProjOp.find(ProjectorType::PROJ)->second.push_back(Matrix(this->mspSetup->fast().at(iL).size(), this->mGrid.size()));
         this->mProjOp.find(ProjectorType::DIVRM0)->second.push_back(Matrix(this->mspSetup->fast().at(iL).size(), this->mGrid.size()));
         this->mProjOp.find(ProjectorType::DIFF)->second.push_back(Matrix(this->mspSetup->fast().at(iL).size(), this->mGrid.size()));
         this->mProjOp.find(ProjectorType::DIFF_PROJ)->second.push_back(Matrix(this->mspSetup->fast().at(iL).size(), this->mGrid.size()));
         this->mProjOp.find(ProjectorType::DIVRDIFFR)->second.push_back(Matrix(this->mspSetup->fast().at(iL).size(), this->mGrid.size()));
         this->mProjOp.find(ProjectorType::LAPLH)->second.push_back(Matrix(this->mspSetup->fast().at(iL).size(), this->mGrid.size()));
         this->mProjOp.find(ProjectorType::LAPLH_DIVRDIFFR)->second.push_back(Matrix(this->mspSetup->fast().at(iL).size(), this->mGrid.size()));
         this->mProjOp.find(ProjectorType::DIVRLAPLHM0)->second.push_back(Matrix(this->mspSetup->fast().at(iL).size(), this->mGrid.size()));
         this->mProjOp.find(ProjectorType::DIFFLAPLH)->second.push_back(Matrix(this->mspSetup->fast().at(iL).size(), this->mGrid.size()));
         this->mProjOp.find(ProjectorType::DIFFLAPLH_DIFFDIVRDIFFR)->second.push_back(Matrix(this->mspSetup->fast().at(iL).size(), this->mGrid.size()));

         op.resize(this->mGrid.size(), this->mspSetup->fast().at(iL).size());

         // Projector: P
         std::map<ProjectorType::Id,std::vector<Matrix> >::iterator projIt = this->mProjOp.find(ProjectorType::PROJ);
         Polynomial::WorlandPolynomial::Wnl(op, ipoly, l, igrid);
         projIt->second.at(iL) = op.transpose();

         // 1/R Projector: 1/R P, set to 0 for m = 0 (always appears combined with \partial_\theta)
         projIt = this->mProjOp.find(ProjectorType::DIVRM0);
         if(l == 0)
         {
            op.setZero();
         } else
         {
            Polynomial::WorlandPolynomial::r_1Wnl(op, itmp, l, igrid);
         }
         projIt->second.at(iL) = op.transpose();

         // Projector: D P
         projIt = this->mProjOp.find(ProjectorType::DIFF);
         Polynomial::WorlandPolynomial::dWnl(op, itmp, l, igrid);
         projIt->second.at(iL) = op.transpose();

         // Projector: D P m > 0, P m = 0
         projIt = this->mProjOp.find(ProjectorType::DIFF_PROJ);
         if(l == 0)
         {
            Polynomial::WorlandPolynomial::Wnl(op, itmp, 1, igrid);
         } else
         {
            Polynomial::WorlandPolynomial::dWnl(op, itmp, l, igrid);
         }
         projIt->second.at(iL) = op.transpose();

         // Projector: 1/R D R P
         projIt = this->mProjOp.find(ProjectorType::DIVRDIFFR);
         Polynomial::WorlandPolynomial::r_1drWnl(op, itmp, l, igrid);
         projIt->second.at(iL) = op.transpose();

         // Projector: horizontal laplacian
         projIt = this->mProjOp.find(ProjectorType::LAPLH);
         Polynomial::WorlandPolynomial::claplhWnl(op, itmp, l, igrid);
         projIt->second.at(iL) = op.transpose();

         // Projector: horizontal laplacian m > 0, 1/r D r m = 0
         projIt = this->mProjOp.find(ProjectorType::LAPLH_DIVRDIFFR);
         if(l == 0)
         {
            Polynomial::WorlandPolynomial::r_1drWnl(op, itmp, 1, igrid);
         } else
         {
            Polynomial::WorlandPolynomial::claplhWnl(op, itmp, l, igrid);
         }
         projIt->second.at(iL) = op.transpose();

         // Projector: 1/R horizontal laplacian, set to 0 for m = 0 (always appears combined with \partial_\theta)
         projIt = this->mProjOp.find(ProjectorType::DIVRLAPLHM0);
         if(l == 0)
         {
            op.setZero();
         } else
         {
            Polynomial::WorlandPolynomial::r_1claplhWnl(op, itmp, l, igrid);
         }
         projIt->second.at(iL) = op.transpose();

         // Projector: diff horizontal laplacian
         projIt = this->mProjOp.find(ProjectorType::DIFFLAPLH);
         Polynomial::WorlandPolynomial::dclaplhWnl(op, itmp, l, igrid);
         projIt->second.at(iL) = op.transpose();

         // Projector: diff horizontal laplacian m > 0, D 1/r D r m =0
         projIt = this->mProjOp.find(ProjectorType::DIFFLAPLH_DIFFDIVRDIFFR);
         if(l == 0)
         {
            Polynomial::WorlandPolynomial::dr_1drWnl(op, itmp, 1, igrid);
         } else
         {
            Polynomial::WorlandPolynomial::dclaplhWnl(op, itmp, l, igrid);
         }
         projIt->second.at(iL) = op.transpose();

         // Set resolution and m
         PyTuple_SetItem(pArgs, 0, PyLong_FromLong(this->mspSetup->fast().at(iL).size()));
         PyTuple_SetItem(pArgs, 1, PyLong_FromLong(l));
         // ... create boundray condition (none)
         pValue = PyDict_New();
         PyDict_SetItem(pValue, PyLong_FromLong(0), PyLong_FromLong(0));
         PyTuple_SetItem(pArgs, 2, pValue);
         // Set coefficient
         PyTuple_SetItem(pArgs, 3, PyFloat_FromDouble(1.0));

         // Call i4
         SparseMatrix matI4(this->mspSetup->fast().at(iL).size(),this->mspSetup->fast().at(iL).size());
         PythonWrapper::setFunction("i4");
         pValue = PythonWrapper::callFunction(pArgs);
         // Fill matrix
         PythonWrapper::fillMatrix(matI4, pValue);
         Py_DECREF(pValue);

         // Call i4divrdiff
         SparseMatrix matI4DIVRDIFF(this->mspSetup->fast().at(iL).size(),this->mspSetup->fast().at(iL).size());
         PythonWrapper::setFunction("i4r_1d1");
         pValue = PythonWrapper::callFunction(pArgs);
         // Fill matrix
         PythonWrapper::fillMatrix(matI4DIVRDIFF, pValue);
         Py_DECREF(pValue);

         // Call i6
         SparseMatrix matI6(this->mspSetup->fast().at(iL).size(),this->mspSetup->fast().at(iL).size());
         PythonWrapper::setFunction("i6");
         pValue = PythonWrapper::callFunction(pArgs);
         // Fill matrix
         PythonWrapper::fillMatrix(matI6, pValue);
         Py_DECREF(pValue);

         // Call i6divrdiff
         SparseMatrix matI6DIVRDIFF(this->mspSetup->fast().at(iL).size(),this->mspSetup->fast().at(iL).size());
         PythonWrapper::setFunction("i6r_1d1");
         pValue = PythonWrapper::callFunction(pArgs);
         // Fill matrix
         PythonWrapper::fillMatrix(matI6DIVRDIFF, pValue);
         Py_DECREF(pValue);

         // Call i6laplh
         SparseMatrix matI6LAPLH(this->mspSetup->fast().at(iL).size(),this->mspSetup->fast().at(iL).size());
         PythonWrapper::setFunction("i6laplh");
         pValue = PythonWrapper::callFunction(pArgs);
         // Fill matrix
         PythonWrapper::fillMatrix(matI6LAPLH, pValue);
         Py_DECREF(pValue);

         SparseMatrix matI2_1;
         SparseMatrix matI4_1;
         SparseMatrix matI4DR;
         if(l == 0)
         {
            // Projection onto radial derivative => l = 1
            PyTuple_SetItem(pArgs, 1, PyLong_FromLong(1));

            // Call i2
            matI2_1.resize(this->mspSetup->fast().at(iL).size(),this->mspSetup->fast().at(iL).size());
            PythonWrapper::setFunction("i2");
            pValue = PythonWrapper::callFunction(pArgs);
            // Fill matrix
            PythonWrapper::fillMatrix(matI2_1, pValue);
            Py_DECREF(pValue);

            // Call i4
            matI4_1.resize(this->mspSetup->fast().at(iL).size(),this->mspSetup->fast().at(iL).size());
            PythonWrapper::setFunction("i4");
            pValue = PythonWrapper::callFunction(pArgs);
            // Fill matrix
            PythonWrapper::fillMatrix(matI4_1, pValue);
            Py_DECREF(pValue);

            // Call i4dr
            matI4DR.resize(this->mspSetup->fast().at(iL).size(),this->mspSetup->fast().at(iL).size());
            PythonWrapper::setFunction("i4dr");
            pValue = PythonWrapper::callFunction(pArgs);
            // Fill matrix
            PythonWrapper::fillMatrix(matI4DR, pValue);
            Py_DECREF(pValue);
         }

         // Allocate memory for the weighted integrator
         this->mIntgOp.find(IntegratorType::INTG)->second.push_back(Matrix(this->mGrid.size(), this->mspSetup->fast().at(iL).size()));
         this->mIntgOp.find(IntegratorType::INTGI4DIVRM0)->second.push_back(Matrix(this->mGrid.size(), this->mspSetup->fast().at(iL).size()));
         this->mIntgOp.find(IntegratorType::INTGI4DIVRDIFFR_I2)->second.push_back(Matrix(this->mGrid.size(), this->mspSetup->fast().at(iL).size()));
         this->mIntgOp.find(IntegratorType::INTGI6DIVRM0)->second.push_back(Matrix(this->mGrid.size(), this->mspSetup->fast().at(iL).size()));
         this->mIntgOp.find(IntegratorType::INTGI6DIVRDIFFR_I4)->second.push_back(Matrix(this->mGrid.size(), this->mspSetup->fast().at(iL).size()));
         this->mIntgOp.find(IntegratorType::INTGI6LAPLH_I4DR)->second.push_back(Matrix(this->mGrid.size(), this->mspSetup->fast().at(iL).size()));

         // Integrator: INTG
         std::map<IntegratorType::Id,std::vector<Matrix> >::iterator intgIt = this->mIntgOp.find(IntegratorType::INTG);
         intgIt->second.at(iL) = (this->mProjOp.find(ProjectorType::PROJ)->second.at(iL)*this->mWeights.asDiagonal()).transpose();
         if(intgIt->second.at(iL).rows() != this->mGrid.size()|| intgIt->second.at(iL).cols() != this->mspSetup->fast().at(iL).size())
         {
            throw Exception("Cylindrical Worland transform operators not setup properly!");
         }

         // Integrator I4DIVRM0
         intgIt =  this->mIntgOp.find(IntegratorType::INTGI4DIVRM0);
         if(l == 0)
         {
            intgIt->second.at(iL).setZero();
         } else
         {
            // Integrator onto W_n^{l-1} basis
            Polynomial::WorlandPolynomial::Wnl(op, ipoly, std::abs(l-1), igrid);
            intgIt->second.at(iL) = (op.transpose()*this->mWeights.asDiagonal()).transpose();
            // Compute 1/r on W_n^{l-1}, integrator onto W_n^{l-1} and apply quasi-inverse
            Polynomial::WorlandPolynomial::r_1Wnl(op, ipoly, std::abs(l-1), igrid);
            // Integrator onto W_n^{l-1} and apply quasi-inverse
            intgIt->second.at(iL) = (matI4*this->mProjOp.find(ProjectorType::PROJ)->second.at(iL)*this->mWeights.asDiagonal()*op*intgIt->second.at(iL).transpose()).transpose();
         }
         if(intgIt->second.at(iL).rows() != this->mGrid.size()|| intgIt->second.at(iL).cols() != this->mspSetup->fast().at(iL).size())
         {
            throw Exception("Cylindrical Worland transform operators not setup properly!");
         }

         // Integrator I4DIVRDIFFR_I2
         intgIt =  this->mIntgOp.find(IntegratorType::INTGI4DIVRDIFFR_I2);
         if(l == 0)
         {
            // Integrator onto W_n^{1} basis
            Polynomial::WorlandPolynomial::Wnl(op, ipoly, 1, igrid);
            intgIt->second.at(iL) = (matI2_1*op.transpose()*this->mWeights.asDiagonal()).transpose();
         } else
         {
            // Mutiply by R in physical space and integrator onto W_n^{l} and apply quasi-inverse
            intgIt->second.at(iL) = (matI4DIVRDIFF*this->mProjOp.find(ProjectorType::PROJ)->second.at(iL)*this->mWeights.asDiagonal()*this->mGrid.asDiagonal()).transpose();
         }
         if(intgIt->second.at(iL).rows() != this->mGrid.size()|| intgIt->second.at(iL).cols() != this->mspSetup->fast().at(iL).size())
         {
            throw Exception("Cylindrical Worland transform operators not setup properly!");
         }

         // Integrator I6DIVRM0
         intgIt =  this->mIntgOp.find(IntegratorType::INTGI6DIVRM0);
         if(l == 0)
         {
            intgIt->second.at(iL).setZero();
         } else
         {
            // Integrator onto W_n^{l-1} basis
            Polynomial::WorlandPolynomial::Wnl(op, ipoly, std::abs(l-1), igrid);
            intgIt->second.at(iL) = (op.transpose()*this->mWeights.asDiagonal()).transpose();
            // Compute 1/r on W_n^{l-1}, integrator onto W_n^{l-1} and apply quasi-inverse
            Polynomial::WorlandPolynomial::r_1Wnl(op, ipoly, std::abs(l-1), igrid);
            // Integrator onto W_n^{l-1} and apply quasi-inverse
            intgIt->second.at(iL) = (matI6*this->mProjOp.find(ProjectorType::PROJ)->second.at(iL)*this->mWeights.asDiagonal()*op*intgIt->second.at(iL).transpose()).transpose();
         }
         if(intgIt->second.at(iL).rows() != this->mGrid.size()|| intgIt->second.at(iL).cols() != this->mspSetup->fast().at(iL).size())
         {
            throw Exception("Cylindrical Worland transform operators not setup properly!");
         }

         // Integrator I6DIVRDIFFR_I4
         intgIt =  this->mIntgOp.find(IntegratorType::INTGI6DIVRDIFFR_I4);
         if(l == 0)
         {
            // Integrator onto W_n^{1} basis
            Polynomial::WorlandPolynomial::Wnl(op, ipoly, 1, igrid);
            intgIt->second.at(iL) = (matI4_1*op.transpose()*this->mWeights.asDiagonal()).transpose();
         } else
         {
            // Multiply by R in physical space and integrator onto W_n^{l} and apply quasi-inverse
            intgIt->second.at(iL) = (matI6DIVRDIFF*this->mProjOp.find(ProjectorType::PROJ)->second.at(iL)*this->mWeights.asDiagonal()*this->mGrid.asDiagonal()).transpose();
         }
         if(intgIt->second.at(iL).rows() != this->mGrid.size()|| intgIt->second.at(iL).cols() != this->mspSetup->fast().at(iL).size())
         {
            throw Exception("Cylindrical Worland transform operators not setup properly!");
         }

         // Integrator I6LAPLH_I4DR
         intgIt =  this->mIntgOp.find(IntegratorType::INTGI6LAPLH_I4DR);
         if(l == 0)
         {
            // Integrator onto W_n^{l} and apply quasi-inverse
            intgIt->second.at(iL) = (matI4DR*this->mProjOp.find(ProjectorType::PROJ)->second.at(iL)*this->mWeights.asDiagonal()).transpose();
         } else
         {
            // Integrator onto W_n^{l} and apply quasi-inverse
            intgIt->second.at(iL) = (matI6LAPLH*this->mProjOp.find(ProjectorType::PROJ)->second.at(iL)*this->mWeights.asDiagonal()).transpose();
         }
         if(intgIt->second.at(iL).rows() != this->mGrid.size()|| intgIt->second.at(iL).cols() != this->mspSetup->fast().at(iL).size())
         {
            throw Exception("Cylindrical Worland transform operators not setup properly!");
         }
      }

      // Cleanup
      PythonWrapper::finalize();
   }

   void CylinderWorlandTransform::integrate(MatrixZ& rSpecVal, const MatrixZ& physVal, CylinderWorlandTransform::IntegratorType::Id integrator)
   {
      // assert right sizes for input matrix
      assert(physVal.rows() == this->mspSetup->fwdSize());
      assert(physVal.cols() == this->mspSetup->howmany());

      // assert right sizes for output matrix
      assert(rSpecVal.cols() == this->mspSetup->howmany());

      // Compute basic integrator
      if(integrator == CylinderWorlandTransform::IntegratorType::INTGINVLAPLH)
      {
         this->setIntegrator(rSpecVal, physVal, this->mIntgOp.find(CylinderWorlandTransform::IntegratorType::INTG)->second);

      } else if(this->mIntgOp.count(integrator) == 1)
      {
         this->setIntegrator(rSpecVal, physVal, this->mIntgOp.find(integrator)->second);

      } else
      {
         throw Exception("Requested an unknown integrator");
      }

   }

   void CylinderWorlandTransform::project(MatrixZ& rPhysVal, const MatrixZ& specVal, CylinderWorlandTransform::ProjectorType::Id projector)
   {
      // assert right sizes for input  matrix
      assert(specVal.cols() == this->mspSetup->howmany());

      // assert right sizes for output matrix
      assert(rPhysVal.rows() == this->mspSetup->fwdSize());
      assert(rPhysVal.cols() == this->mspSetup->howmany());

      // Compute projection
      if(this->mProjOp.count(projector) == 1)
      {
         this->setProjector(rPhysVal, specVal, this->mProjOp.find(projector)->second);

      } else
      {
         throw Exception("Requested an unknown projector");
      }
   }

   void CylinderWorlandTransform::setIntegrator(MatrixZ& rSpecVal, const MatrixZ& physVal, const std::vector<Matrix>& ops)
   {
      // Compute integration
      int start = 0;
      int physRows = this->mspSetup->fwdSize(); 
      for(size_t i = 0; i < ops.size(); i++)
      {
         int cols = this->mspSetup->mult()(i);
         int specRows = ops.at(i).cols();
         rSpecVal.block(0, start, specRows, cols) = ops.at(i).transpose()*physVal.block(0,start, physRows, cols);
         start += cols;
      }
   }

   void CylinderWorlandTransform::setProjector(MatrixZ& rPhysVal, const MatrixZ& specVal, const std::vector<Matrix>& ops)
   {
      int start = 0;
      int physRows = this->mspSetup->fwdSize(); 
      for(size_t i = 0; i < ops.size(); i++)
      {
         int cols = this->mspSetup->mult()(i);
         int specRows = ops.at(i).rows();
         rPhysVal.block(0, start, physRows, cols) = ops.at(i).transpose()*specVal.block(0,start, specRows, cols);
         start += cols;
      }
   }

#ifdef QUICC_STORAGEPROFILE
   MHDFloat CylinderWorlandTransform::requiredStorage() const
   {
      MHDFloat mem = 0.0;

      // Storage for the grid and weight
      mem += static_cast<MHDFloat>(Debug::MemorySize<MHDFloat>::BYTES)*this->mGrid.size();
      mem += static_cast<MHDFloat>(Debug::MemorySize<MHDFloat>::BYTES)*this->mWeights.size();

      // Storage for the projectors
      for(std::map<ProjectorType::Id,std::vector<Matrix> >::const_iterator projIt = this->mProjOp.begin(); projIt != this->mProjOp.end(); ++projIt)
      {
         for(size_t i = 0; i < projIt->second.size(); i++)
         {
            mem += static_cast<MHDFloat>(Debug::MemorySize<MHDFloat>::BYTES)*projIt->second.at(i).size();
         }
      }

      // Storage for the projectors
      for(std::map<IntegratorType::Id,std::vector<Matrix> >::const_iterator intgIt = this->mIntgOp.begin(); intgIt != this->mIntgOp.end(); ++intgIt)
      {
         for(size_t i = 0; i < intgIt->second.size(); i++)
         {
            mem += static_cast<MHDFloat>(Debug::MemorySize<MHDFloat>::BYTES)*intgIt->second.at(i).size();
         }
      }

      return mem;
   }
#endif // QUICC_STORAGEPROFILE

}
}
