/** 
 * @file SphereWorlandTransform.cpp
 * @brief Source of the implementation of the Worland transform in a sphere
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "PolynomialTransforms/SphereWorlandTransform.hpp"

// Project includes
//
#include "Exceptions/Exception.hpp"
#include "Base/MathConstants.hpp"
#include "Quadratures/WorlandChebyshevRule.hpp"
#include "Quadratures/LegendreRule.hpp"
#include "PolynomialTransforms/WorlandPolynomial.hpp"
#include "Python/PythonWrapper.hpp"

namespace GeoMHDiSCC {

namespace Transform {

   Array SphereWorlandTransform::generateGrid(const int size)
   {
      // Initialise grid storage
      Array grid(size);
      Array tmp(size);

      WorlandChebyshevRule::computeQuadrature(grid, tmp, size);

      return grid;
   }

   SphereWorlandTransform::SphereWorlandTransform()
   {
      // Initialise the Python interpreter wrapper
      PythonWrapper::init();
   }

   SphereWorlandTransform::~SphereWorlandTransform()
   {
   }

   void SphereWorlandTransform::init(SphereWorlandTransform::SharedSetupType spSetup)
   {
      // Store the shared pointer to setup object
      this->mspSetup = spSetup;

      // Initialise the quadrature grid and weights and operators
      this->initOperators();
   }

   void SphereWorlandTransform::requiredOptions(std::set<NonDimensional::Id>& list, const Dimensions::Transform::Id dimId) const
   {
      //
      // No possible options
      //
   }

   void SphereWorlandTransform::setOptions(const std::map<NonDimensional::Id, MHDFloat>& options, const Dimensions::Transform::Id dimId)
   {
      //
      // No possible options
      //
   }

   const Array& SphereWorlandTransform::meshGrid() const
   {
      if(this->mGrid.size() == 0 || this->mWeights.size() == 0)
      {
         throw Exception("SphereWorlandTransform has not been initialised!");
      }

      return this->mGrid;
   }

   void SphereWorlandTransform::initOperators()
   {
      // Initialise python wrapper
      PythonWrapper::import("quicc.geometry.spherical.sphere_radius_worland");

      this->mGrid.resize(this->mspSetup->fwdSize());
      this->mWeights.resize(this->mspSetup->fwdSize());

      // Set the grid and weights
      internal::Array igrid, iweights;
      WorlandChebyshevRule::computeQuadrature(this->mGrid, this->mWeights, igrid, iweights, this->mspSetup->fwdSize());

      // Create grids for energy calculations
      Array legGrid(std::ceil(4.0*this->mspSetup->fwdSize()/3.) + 1);
      Array legWeights(legGrid.size());
      internal::Array ilegGrid, ilegWeights;
      LegendreRule::computeQuadrature(legGrid, legWeights, ilegGrid, ilegWeights, legGrid.size());
      legGrid.array() = ((legGrid.array() + 1.0)/2.0);
      ilegGrid.array() = ((ilegGrid.array() + 1.0)/2.0);

      Array nrgGrid(std::ceil(4.0*this->mspSetup->fwdSize()/3.) + 1);
      Array nrgWeights(nrgGrid.size());
      internal::Array inrgGrid, inrgWeights;
      WorlandChebyshevRule::computeQuadrature(nrgGrid, nrgWeights, inrgGrid, inrgWeights, std::ceil(4.0*this->mspSetup->fwdSize()/3.) + 1);

      // Reserve storage for the projectors, 1st derivative
      this->mProjOp.insert(std::make_pair(ProjectorType::PROJ,std::vector<Matrix>()));
      this->mProjOp.find(ProjectorType::PROJ)->second.reserve(this->mspSetup->slow().size());
      this->mProjOp.insert(std::make_pair(ProjectorType::DIVR,std::vector<Matrix>()));
      this->mProjOp.find(ProjectorType::DIVR)->second.reserve(this->mspSetup->slow().size());
      this->mProjOp.insert(std::make_pair(ProjectorType::DIFF,std::vector<Matrix>()));
      this->mProjOp.find(ProjectorType::DIFF)->second.reserve(this->mspSetup->slow().size());
      this->mProjOp.insert(std::make_pair(ProjectorType::DIFFR,std::vector<Matrix>()));
      this->mProjOp.find(ProjectorType::DIFFR)->second.reserve(this->mspSetup->slow().size());
      this->mProjOp.insert(std::make_pair(ProjectorType::DIVRDIFFR,std::vector<Matrix>()));
      this->mProjOp.find(ProjectorType::DIVRDIFFR)->second.reserve(this->mspSetup->slow().size());
      this->mProjOp.insert(std::make_pair(ProjectorType::SLAPL,std::vector<Matrix>()));
      this->mProjOp.find(ProjectorType::SLAPL)->second.reserve(this->mspSetup->slow().size());
      this->mProjOp.insert(std::make_pair(ProjectorType::ENERGY_PROJ,std::vector<Matrix>()));
      this->mProjOp.find(ProjectorType::ENERGY_PROJ)->second.reserve(this->mspSetup->slow().size());
      this->mProjOp.insert(std::make_pair(ProjectorType::ENERGY_DIFFR,std::vector<Matrix>()));
      this->mProjOp.find(ProjectorType::ENERGY_DIFFR)->second.reserve(this->mspSetup->slow().size());

      // Reserve storage for the weighted projectors 
      this->mIntgOp.insert(std::make_pair(IntegratorType::INTG,std::vector<Matrix>()));
      this->mIntgOp.find(IntegratorType::INTG)->second.reserve(this->mspSetup->slow().size());
      this->mIntgOp.insert(std::make_pair(IntegratorType::INTGR,std::vector<Matrix>()));
      this->mIntgOp.find(IntegratorType::INTGR)->second.reserve(this->mspSetup->slow().size());
      this->mIntgOp.insert(std::make_pair(IntegratorType::INTGQ4,std::vector<Matrix>()));
      this->mIntgOp.find(IntegratorType::INTGQ4)->second.reserve(this->mspSetup->slow().size());
      this->mIntgOp.insert(std::make_pair(IntegratorType::INTGS4,std::vector<Matrix>()));
      this->mIntgOp.find(IntegratorType::INTGS4)->second.reserve(this->mspSetup->slow().size());
      this->mIntgOp.insert(std::make_pair(IntegratorType::INTGT,std::vector<Matrix>()));
      this->mIntgOp.find(IntegratorType::INTGT)->second.reserve(this->mspSetup->slow().size());
      this->mIntgOp.insert(std::make_pair(IntegratorType::INTGQ2,std::vector<Matrix>()));
      this->mIntgOp.find(IntegratorType::INTGQ2)->second.reserve(this->mspSetup->slow().size());
      this->mIntgOp.insert(std::make_pair(IntegratorType::INTGS2,std::vector<Matrix>()));
      this->mIntgOp.find(IntegratorType::INTGS2)->second.reserve(this->mspSetup->slow().size());
      this->mIntgOp.insert(std::make_pair(IntegratorType::ENERGY_INTG,std::vector<Matrix>()));
      this->mIntgOp.find(IntegratorType::ENERGY_INTG)->second.reserve(this->mspSetup->slow().size());
      this->mIntgOp.insert(std::make_pair(IntegratorType::ENERGY_R2,std::vector<Matrix>()));
      this->mIntgOp.find(IntegratorType::ENERGY_R2)->second.reserve(this->mspSetup->slow().size());

      // Prepare arguments to Python matrices call
      PyObject *pArgs, *pValue;
      pArgs = PyTuple_New(4);

      // Loop over harmonic degrees
      Matrix op;
      internal::Matrix  itmp, ipoly;
      for(int iL = 0; iL < this->mspSetup->slow().size(); iL++)
      {
         int l = this->mspSetup->slow()(iL);

         op.resize(this->mGrid.size(), this->mspSetup->fast().at(iL).size());

         // Projector: P
         std::map<ProjectorType::Id,std::vector<Matrix> >::iterator projIt = this->mProjOp.find(ProjectorType::PROJ);
         projIt->second.push_back(Matrix(this->mspSetup->fast().at(iL).size(), this->mGrid.size()));
         Polynomial::WorlandPolynomial::Wnl(op, ipoly, l, igrid);
         projIt->second.at(iL) = op.transpose();

         // 1/R Projector: 1/R P
         projIt = this->mProjOp.find(ProjectorType::DIVR);
         projIt->second.push_back(Matrix(this->mspSetup->fast().at(iL).size(), this->mGrid.size()));
         Polynomial::WorlandPolynomial::r_1Wnl(op, itmp, l, igrid);
         projIt->second.at(iL) = op.transpose();

         // First derivative: D
         projIt = this->mProjOp.find(ProjectorType::DIFF);
         projIt->second.push_back(Matrix(this->mspSetup->fast().at(iL).size(), this->mGrid.size()));
         Polynomial::WorlandPolynomial::dWnl(op, itmp, l, igrid);
         projIt->second.at(iL) = op.transpose();

         // D R
         projIt = this->mProjOp.find(ProjectorType::DIFFR);
         projIt->second.push_back(Matrix(this->mspSetup->fast().at(iL).size(), this->mGrid.size()));
         Polynomial::WorlandPolynomial::drWnl(op, itmp, l, igrid);
         projIt->second.at(iL) = op.transpose();

         // 1/R D R
         projIt = this->mProjOp.find(ProjectorType::DIVRDIFFR);
         projIt->second.push_back(Matrix(this->mspSetup->fast().at(iL).size(), this->mGrid.size()));
         Polynomial::WorlandPolynomial::r_1drWnl(op, itmp, l, igrid);
         projIt->second.at(iL) = op.transpose();

         // Spherical laplacian: D^2 + 2/R D - l(l+1)/R^2
         projIt = this->mProjOp.find(ProjectorType::SLAPL);
         projIt->second.push_back(Matrix(this->mspSetup->fast().at(iL).size(), this->mGrid.size()));
         Polynomial::WorlandPolynomial::slaplWnl(op, itmp, l, igrid);
         projIt->second.at(iL) = op.transpose();

         // Set resolution and l
         PyTuple_SetItem(pArgs, 0, PyLong_FromLong(this->mspSetup->fast().at(iL).size()));
         PyTuple_SetItem(pArgs, 1, PyLong_FromLong(l));
         // ... create boundray condition (none)
         pValue = PyDict_New();
         PyDict_SetItem(pValue, PyLong_FromLong(0), PyLong_FromLong(0));
         PyTuple_SetItem(pArgs, 2, pValue);
         // Set coefficient
         PyTuple_SetItem(pArgs, 3, PyFloat_FromDouble(1.0));

         // Call i2
         SparseMatrix matI2(this->mspSetup->fast().at(iL).size(),this->mspSetup->fast().at(iL).size());
         PythonWrapper::setFunction("i2");
         pValue = PythonWrapper::callFunction(pArgs);
         // Fill matrix
         PythonWrapper::fillMatrix(matI2, pValue);
         Py_DECREF(pValue);

         // Call i4
         SparseMatrix matI4(this->mspSetup->fast().at(iL).size(),this->mspSetup->fast().at(iL).size());
         PythonWrapper::setFunction("i4");
         pValue = PythonWrapper::callFunction(pArgs);
         // Fill matrix
         PythonWrapper::fillMatrix(matI4, pValue);
         Py_DECREF(pValue);

         // Call r2
         int energySize = 2*this->mspSetup->fast().at(iL).size();
         PyTuple_SetItem(pArgs, 0, PyLong_FromLong(energySize+1));
         PyTuple_SetItem(pArgs, 1, PyLong_FromLong(2*l));
         SparseMatrix matTmp(energySize+1,energySize+1);
         PythonWrapper::setFunction("r2");
         pValue = PythonWrapper::callFunction(pArgs);
         // Fill matrix
         PythonWrapper::fillMatrix(matTmp, pValue);
         SparseMatrix matNRG_R2(energySize+1,energySize);
         matNRG_R2 = matTmp.leftCols(energySize);
         matTmp.resize(0,0);
         Py_DECREF(pValue);

         // Integrator: INTG
         std::map<IntegratorType::Id,std::vector<Matrix> >::iterator intgIt = this->mIntgOp.find(IntegratorType::INTG);
         intgIt->second.push_back(Matrix(this->mGrid.size(), this->mspSetup->fast().at(iL).size()));
         intgIt->second.at(iL) = (this->mProjOp.find(ProjectorType::PROJ)->second.at(iL)*this->mWeights.asDiagonal()).transpose();
         if(intgIt->second.at(iL).rows() != this->mGrid.size()|| intgIt->second.at(iL).cols() != this->mspSetup->fast().at(iL).size())
         {
            throw Exception("Spherical Worland transform operators not setup properly!");
         }

         // Integrator: INTGR
         intgIt = this->mIntgOp.find(IntegratorType::INTGR);
         intgIt->second.push_back(Matrix(this->mGrid.size(), this->mspSetup->fast().at(iL).size()));
         intgIt->second.at(iL) = (this->mProjOp.find(ProjectorType::PROJ)->second.at(iL)*this->mWeights.asDiagonal()*this->mGrid.asDiagonal()).transpose();
         if(intgIt->second.at(iL).rows() != this->mGrid.size()|| intgIt->second.at(iL).cols() != this->mspSetup->fast().at(iL).size())
         {
            throw Exception("Spherical Worland transform operators not setup properly!");
         }

         // Integrator: INTGQ4
         intgIt = this->mIntgOp.find(IntegratorType::INTGQ4);
         intgIt->second.push_back(Matrix(this->mGrid.size(), this->mspSetup->fast().at(iL).size()));
         if(l == 0)
         {
            intgIt->second.at(iL).setZero();
         } else
         {
            Polynomial::WorlandPolynomial::Wnl(op, ipoly, std::abs(l-1), igrid);
            intgIt->second.at(iL) = (op.transpose()*this->mWeights.asDiagonal()).transpose();
            Polynomial::WorlandPolynomial::r_1Wnl(op, ipoly, std::abs(l-1), igrid);
            intgIt->second.at(iL) = (matI4*this->mProjOp.find(ProjectorType::PROJ)->second.at(iL)*this->mWeights.asDiagonal()*op*intgIt->second.at(iL).transpose()).transpose();
         }
         if(intgIt->second.at(iL).rows() != this->mGrid.size()|| intgIt->second.at(iL).cols() != this->mspSetup->fast().at(iL).size())
         {
            throw Exception("Spherical Worland transform operators not setup properly!");
         }

         // Integrator: INTGS4
         intgIt = this->mIntgOp.find(IntegratorType::INTGS4);
         intgIt->second.push_back(Matrix(this->mGrid.size(), this->mspSetup->fast().at(iL).size()));
         if(l == 0)
         {
            intgIt->second.at(iL).setZero();
         } else
         {
            Polynomial::WorlandPolynomial::Wnl(op, ipoly, std::abs(l-1), igrid);
            intgIt->second.at(iL) = (op.transpose()*this->mWeights.asDiagonal()).transpose();
            Polynomial::WorlandPolynomial::r_1drWnl(op, ipoly, std::abs(l-1), igrid);
            intgIt->second.at(iL) = (matI4*this->mProjOp.find(ProjectorType::PROJ)->second.at(iL)*this->mWeights.asDiagonal()*op*intgIt->second.at(iL).transpose()).transpose();
         }
         if(intgIt->second.at(iL).rows() != this->mGrid.size()|| intgIt->second.at(iL).cols() != this->mspSetup->fast().at(iL).size())
         {
            throw Exception("Spherical Worland transform operators not setup properly!");
         }

         // Integrator: INTGT
         intgIt = this->mIntgOp.find(IntegratorType::INTGT);
         intgIt->second.push_back(Matrix(this->mGrid.size(), this->mspSetup->fast().at(iL).size()));
         if(l == 0)
         {
            intgIt->second.at(iL).setZero();
         } else
         {
            intgIt->second.at(iL) = (matI2*this->mProjOp.find(ProjectorType::PROJ)->second.at(iL)*this->mWeights.asDiagonal()).transpose();
         }
         if(intgIt->second.at(iL).rows() != this->mGrid.size()|| intgIt->second.at(iL).cols() != this->mspSetup->fast().at(iL).size())
         {
            throw Exception("Spherical Worland transform operators not setup properly!");
         }

         // Integrator: INTGQ2
         intgIt = this->mIntgOp.find(IntegratorType::INTGQ2);
         intgIt->second.push_back(Matrix(this->mGrid.size(), this->mspSetup->fast().at(iL).size()));
         if(l == 0)
         {
            intgIt->second.at(iL).setZero();
         } else
         {
            Polynomial::WorlandPolynomial::Wnl(op, ipoly, std::abs(l-1), igrid);
            intgIt->second.at(iL) = (op.transpose()*this->mWeights.asDiagonal()).transpose();
            Polynomial::WorlandPolynomial::r_1Wnl(op, ipoly, std::abs(l-1), igrid);
            intgIt->second.at(iL) = (matI2*this->mProjOp.find(ProjectorType::PROJ)->second.at(iL)*this->mWeights.asDiagonal()*op*intgIt->second.at(iL).transpose()).transpose();
         }
         if(intgIt->second.at(iL).rows() != this->mGrid.size()|| intgIt->second.at(iL).cols() != this->mspSetup->fast().at(iL).size())
         {
            throw Exception("Spherical Worland transform operators not setup properly!");
         }

         // Integrator: INTGS2
         intgIt = this->mIntgOp.find(IntegratorType::INTGS2);
         intgIt->second.push_back(Matrix(this->mGrid.size(), this->mspSetup->fast().at(iL).size()));
         if(l == 0)
         {
            intgIt->second.at(iL).setZero();
         } else
         {
            Polynomial::WorlandPolynomial::Wnl(op, ipoly, std::abs(l-1), igrid);
            intgIt->second.at(iL) = (op.transpose()*this->mWeights.asDiagonal()).transpose();
            Polynomial::WorlandPolynomial::r_1drWnl(op, ipoly, std::abs(l-1), igrid);
            intgIt->second.at(iL) = (matI2*this->mProjOp.find(ProjectorType::PROJ)->second.at(iL)*this->mWeights.asDiagonal()*op*intgIt->second.at(iL).transpose()).transpose();
         }
         if(intgIt->second.at(iL).rows() != this->mGrid.size()|| intgIt->second.at(iL).cols() != this->mspSetup->fast().at(iL).size())
         {
            throw Exception("Spherical Worland transform operators not setup properly!");
         }

         //
         // Energy operators
         //
         op.resize(legGrid.size(), energySize+1);
         Polynomial::WorlandPolynomial::Wnl(op, ipoly, 2*l, ilegGrid);
         Array energyWeights = 0.5*(legWeights.asDiagonal()*op).colwise().sum().transpose();

         // Energy projector's size
         op.resize(nrgGrid.size(), this->mspSetup->fast().at(iL).size());

         // Projector
         projIt = this->mProjOp.find(ProjectorType::ENERGY_PROJ);
         projIt->second.push_back(Matrix(this->mspSetup->fast().at(iL).size(), nrgGrid.size()));
         Polynomial::WorlandPolynomial::Wnl(op, itmp, l, inrgGrid);
         projIt->second.at(iL) = op.transpose();

         // Projector: D R
         projIt = this->mProjOp.find(ProjectorType::ENERGY_DIFFR);
         projIt->second.push_back(Matrix(this->mspSetup->fast().at(iL).size(), nrgGrid.size()));
         Polynomial::WorlandPolynomial::drWnl(op, itmp, l, inrgGrid);
         projIt->second.at(iL) = op.transpose();

         // Energy integrator's size
         op.resize(nrgGrid.size(), energySize);

         // Integrator: ENERGY_INTG
         intgIt = this->mIntgOp.find(IntegratorType::ENERGY_INTG);
         intgIt->second.push_back(Matrix(nrgGrid.size(), 1));
         op.resize(nrgGrid.size(), energySize);
         Polynomial::WorlandPolynomial::Wnl(op, ipoly, 2*l, inrgGrid);
         intgIt->second.at(iL).col(0) = (energyWeights.topRows(energySize).transpose()*op.transpose()*nrgWeights.asDiagonal()).transpose();
         if(intgIt->second.at(iL).rows() != nrgGrid.size()|| intgIt->second.at(iL).cols() != 1)
         {
            throw Exception("Spherical Worland transform operators not setup properly!");
         }

         // Integrator: ENERGY_R2
         intgIt = this->mIntgOp.find(IntegratorType::ENERGY_R2);
         intgIt->second.push_back(Matrix(nrgGrid.size(), 1));
         op.resize(nrgGrid.size(), energySize);
         Polynomial::WorlandPolynomial::Wnl(op, ipoly, 2*l, inrgGrid);
         intgIt->second.at(iL).col(0) = (energyWeights.transpose()*matNRG_R2*op.transpose()*nrgWeights.asDiagonal()).transpose();
         if(intgIt->second.at(iL).rows() != nrgGrid.size()|| intgIt->second.at(iL).cols() != 1)
         {
            throw Exception("Spherical Worland transform operators not setup properly!");
         }
      }
   }

   void SphereWorlandTransform::integrate(MatrixZ& rSpecVal, const MatrixZ& physVal, SphereWorlandTransform::IntegratorType::Id integrator)
   {
      // assert right sizes for input matrix
      assert(physVal.rows() == this->mspSetup->fwdSize());
      assert(physVal.cols() == this->mspSetup->howmany());

      // assert right sizes for output matrix
      assert(rSpecVal.cols() == this->mspSetup->howmany());

      // Compute integration
      if(this->mIntgOp.count(integrator) == 1)
      {
         this->setIntegrator(rSpecVal, physVal, this->mIntgOp.find(integrator)->second);

      } else
      {
         throw Exception("Requested an unknown integrator");
      }

   }

   void SphereWorlandTransform::integrate_energy(Array& spectrum, const MatrixZ& specVal, SphereWorlandTransform::ProjectorType::Id projector, SphereWorlandTransform::IntegratorType::Id integrator)
   {
      // assert right sizes for output matrix
      assert(specVal.cols() == this->mspSetup->howmany());
      spectrum.resize(this->mspSetup->howmany());

      // Compute energy integration
      if((projector == SphereWorlandTransform::ProjectorType::ENERGY_PROJ || projector == SphereWorlandTransform::ProjectorType::ENERGY_DIFFR) && (integrator == SphereWorlandTransform::IntegratorType::ENERGY_INTG || integrator == SphereWorlandTransform::IntegratorType::ENERGY_R2))
      {
         this->setEnergyIntegrator(spectrum, specVal, this->mProjOp.find(projector)->second, this->mIntgOp.find(integrator)->second);

      } else
      {
         throw Exception("Requested an unknown energy integrator");
      }

   }

   void SphereWorlandTransform::project(MatrixZ& rPhysVal, const MatrixZ& specVal, SphereWorlandTransform::ProjectorType::Id projector)
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

   void SphereWorlandTransform::setIntegrator(MatrixZ& rSpecVal, const MatrixZ& physVal, const std::vector<Matrix>& ops)
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

   void SphereWorlandTransform::setProjector(MatrixZ& rPhysVal, const MatrixZ& specVal, const std::vector<Matrix>& ops)
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

   void SphereWorlandTransform::setEnergyIntegrator(Array& spectrum, const MatrixZ& specVal, const std::vector<Matrix>& projOps, const std::vector<Matrix>& intgOps)
   {
      // Compute integration
      int start = 0;
      for(size_t i = 0; i < projOps.size(); i++)
      {
         int cols = this->mspSetup->mult()(i);
         int specRows = projOps.at(i).rows();
         MatrixZ tmp = projOps.at(i).transpose()*specVal.block(0,start, specRows, cols);
         tmp.array() = tmp.array()*tmp.conjugate().array();
         spectrum.segment(start, cols).transpose() = intgOps.at(i).transpose()*tmp.real();
         start += cols;
      }
   }

#ifdef QUICC_STORAGEPROFILE
   MHDFloat SphereWorlandTransform::requiredStorage() const
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
