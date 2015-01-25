/** 
 * @file AssociatedLegendreTransform.cpp
 * @brief Source of the implementation of the associated Legendre transform
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// System includes
//

// External includes
//

// Class include
//
#include "PolynomialTransforms/AssociatedLegendreTransform.hpp"

// Project includes
//
#include "Exceptions/Exception.hpp"
#include "Base/MathConstants.hpp"
#include "Quadratures/LegendreRule.hpp"
#include "PolynomialTransforms/AssociatedLegendrePolynomial.hpp"

#include <iostream>
namespace GeoMHDiSCC {

namespace Transform {

   Array AssociatedLegendreTransform::generateGrid(const int size)
   {
      // Initialise grid storage
      Array grid(size);
      Array tmp(size);

      LegendreRule::computeQuadrature(grid, tmp, size);

      grid = grid.array().acos();

      return grid;
   }

   AssociatedLegendreTransform::AssociatedLegendreTransform()
   {
   }

   AssociatedLegendreTransform::~AssociatedLegendreTransform()
   {
   }

   void AssociatedLegendreTransform::init(AssociatedLegendreTransform::SharedSetupType spSetup)
   {
      // Store the shared pointer to setup object
      this->mspSetup = spSetup;

      // Initialise the quadrature grid and weights
      this->initOperators();
   }

   void AssociatedLegendreTransform::requiredOptions(std::set<NonDimensional::Id>& list, const Dimensions::Transform::Id dimId) const
   {
      //
      // No possible options
      //
   }

   void AssociatedLegendreTransform::setOptions(const std::map<NonDimensional::Id, MHDFloat>& options, const Dimensions::Transform::Id dimId)
   {
      //
      // No possible options
      //
   }

   const Array& AssociatedLegendreTransform::meshGrid() const
   {
      if(this->mXGrid.size() == 0 || this->mThGrid.size() == 0 || this->mWeights.size() == 0)
      {
         throw Exception("AssociatedLegendreTransform has not been initialised!");
      }

      return this->mThGrid;
   }

   void AssociatedLegendreTransform::initOperators()
   {
      this->mXGrid.resize(this->mspSetup->fwdSize());
      this->mThGrid.resize(this->mspSetup->fwdSize());
      this->mWeights.resize(this->mspSetup->fwdSize());

      internal::Array igrid, iweights;

      LegendreRule::computeQuadrature(this->mXGrid, this->mWeights, igrid, iweights, this->mspSetup->fwdSize());

      this->mThGrid = this->mXGrid.array().acos();

      // Normalise weights by 2*pi for spherical harmonics
      this->mWeights.array() *= 2.0*Math::PI;

      // Initialise storage for l(l+1) factor
      this->mLl = Array::LinSpaced(this->mspSetup->specSize(), 0, this->mspSetup->specSize()-1);
      this->mLl = this->mLl.array()*(this->mLl.array() + 1.0);
      this->mLl2 = this->mLl.array().pow(2);
      this->mDivLl = this->mLl.array().pow(-1);
      this->mDivLl(0) = 0.0;

      // Reserve storage for the projectors, 1/sin projectors and derivative
      this->mProj.reserve(this->mspSetup->slow().size());
      this->mDivSin.reserve(this->mspSetup->slow().size());
      this->mDiff.reserve(this->mspSetup->slow().size());

      // Loop over harmonic orders
      Matrix op;
      internal::Matrix  itmp, ipoly;
      for(int iM = 0; iM < this->mspSetup->slow().size(); iM++)
      {
         int m = this->mspSetup->slow()(iM);
         int maxL = this->mspSetup->fast().at(iM)(this->mspSetup->fast().at(iM).size()-1);

         // Allocate memory for the projector, 1/sin projector and derivatives
         this->mProj.push_back(Matrix(this->mThGrid.size(), this->mspSetup->fast().at(iM).size()));
         this->mDivSin.push_back(Matrix(this->mThGrid.size(), this->mspSetup->fast().at(iM).size()));
         this->mDiff.push_back(Matrix(this->mThGrid.size(), this->mspSetup->fast().at(iM).size()));

         op.resize(this->mThGrid.size(), maxL - m + 1);

         // Loop over harmonic degrees for projectors
         Polynomial::AssociatedLegendrePolynomial::Plm(op, ipoly, m, igrid);
         for(int iL = 0; iL < this->mspSetup->fast().at(iM).size(); iL++)
         {
            int l = this->mspSetup->fast().at(iM)(iL);

            this->mProj.at(iM).col(iL) = op.col(l - m);
         }

         // Loop over harmonic degrees for derivative
         Polynomial::AssociatedLegendrePolynomial::dPlm(op, itmp, m, ipoly, igrid);
         for(int iL = 0; iL < this->mspSetup->fast().at(iM).size(); iL++)
         {
            int l = this->mspSetup->fast().at(iM)(iL);

            this->mDiff.at(iM).col(iL) = op.col(l - m);
         }

         // Loop over harmonic degrees for 1/\sin\theta
         Polynomial::AssociatedLegendrePolynomial::sin_1Plm(op, itmp, m, igrid);
         for(int iL = 0; iL < this->mspSetup->fast().at(iM).size(); iL++)
         {
            int l = this->mspSetup->fast().at(iM)(iL);

            this->mDivSin.at(iM).col(iL) = op.col(l - m);
         }
      }

   }

   void AssociatedLegendreTransform::integrate(MatrixZ& rSpecVal, const MatrixZ& physVal, AssociatedLegendreTransform::IntegratorType::Id integrator, Arithmetics::Id arithId)
   {
      assert(arithId == Arithmetics::SET);

      // assert right sizes for input matrix
      assert(physVal.rows() == this->mspSetup->fwdSize());
      assert(physVal.cols() == this->mspSetup->howmany());

      // assert right sizes for output matrix
      assert(rSpecVal.cols() == this->mspSetup->howmany());

      // Compute first derivative integration
      if(integrator == AssociatedLegendreTransform::IntegratorType::INTGDIFF)
      {
         // Compute integration of derivative
         int start = 0;
         int physRows = this->mspSetup->fwdSize(); 
         for(size_t i = 0; i < this->mDiff.size(); i++)
         {
            int cols = this->mspSetup->mult()(i);
            int specRows = this->mDiff.at(i).cols();
            rSpecVal.block(0, start, specRows, cols) = this->mDiff.at(i).transpose()*this->mWeights.asDiagonal()*physVal.block(0,start, physRows, cols);
            start += cols;
         }

      } else if(integrator == AssociatedLegendreTransform::IntegratorType::INTGLLDIFF)
      { 
         // Compute integration of derivative
         int start = 0;
         int physRows = this->mspSetup->fwdSize(); 
         for(size_t i = 0; i < this->mDiff.size(); i++)
         {
            int cols = this->mspSetup->mult()(i);
            int specRows = this->mDiff.at(i).cols();
            rSpecVal.block(0, start, specRows, cols) = this->mLl.bottomRows(specRows).asDiagonal()*this->mDiff.at(i).transpose()*this->mWeights.asDiagonal()*physVal.block(0,start, physRows, cols);
            start += cols;
         }

      } else if(integrator == AssociatedLegendreTransform::IntegratorType::INTGDIVLLDIFF)
      { 
         // Compute integration of derivative
         int start = 0;
         int physRows = this->mspSetup->fwdSize(); 
         for(size_t i = 0; i < this->mDiff.size(); i++)
         {
            int cols = this->mspSetup->mult()(i);
            int specRows = this->mDiff.at(i).cols();
            rSpecVal.block(0, start, specRows, cols) = this->mDivLl.bottomRows(specRows).asDiagonal()*this->mDiff.at(i).transpose()*this->mWeights.asDiagonal()*physVal.block(0,start, physRows, cols);
            start += cols;
         }

      } else if(integrator == AssociatedLegendreTransform::IntegratorType::INTGLLDIVSIN)
      { 
         // Compute integration of 1/sin
         int start = 0;
         int physRows = this->mspSetup->fwdSize(); 
         for(size_t i = 0; i < this->mDivSin.size(); i++)
         {
            int cols = this->mspSetup->mult()(i);
            int specRows = this->mDivSin.at(i).cols();
            rSpecVal.block(0, start, specRows, cols) = this->mLl.bottomRows(specRows).asDiagonal()*this->mDivSin.at(i).transpose()*this->mWeights.asDiagonal()*physVal.block(0,start, physRows, cols);
            start += cols;
         }

      } else if(integrator == AssociatedLegendreTransform::IntegratorType::INTGDIVLLDIVSIN)
      { 
         // Compute integration of 1/sin
         int start = 0;
         int physRows = this->mspSetup->fwdSize(); 
         for(size_t i = 0; i < this->mDivSin.size(); i++)
         {
            int cols = this->mspSetup->mult()(i);
            int specRows = this->mDivSin.at(i).cols();
            rSpecVal.block(0, start, specRows, cols) = this->mDivLl.bottomRows(specRows).asDiagonal()*this->mDivSin.at(i).transpose()*this->mWeights.asDiagonal()*physVal.block(0,start, physRows, cols);
            start += cols;
         }

      } else if(integrator == AssociatedLegendreTransform::IntegratorType::INTGDIVSIN)
      { 
         // Compute integration of 1/sin
         int start = 0;
         int physRows = this->mspSetup->fwdSize(); 
         for(size_t i = 0; i < this->mDivSin.size(); i++)
         {
            int cols = this->mspSetup->mult()(i);
            int specRows = this->mDivSin.at(i).cols();
            rSpecVal.block(0, start, specRows, cols) = this->mDivSin.at(i).transpose()*this->mWeights.asDiagonal()*physVal.block(0,start, physRows, cols);
            start += cols;
         }

      } else if(integrator == AssociatedLegendreTransform::IntegratorType::INTGLL2)
      { 
         // Compute integration
         int start = 0;
         int physRows = this->mspSetup->fwdSize(); 
         for(size_t i = 0; i < this->mDiff.size(); i++)
         {
            int cols = this->mspSetup->mult()(i);
            int specRows = this->mProj.at(i).cols();
            rSpecVal.block(0, start, specRows, cols) = this->mLl2.bottomRows(specRows).asDiagonal()*this->mProj.at(i).transpose()*this->mWeights.asDiagonal()*physVal.block(0,start, physRows, cols);
            start += cols;
         }

      } else if(integrator == AssociatedLegendreTransform::IntegratorType::INTGLL)
      { 
         // Compute integration
         int start = 0;
         int physRows = this->mspSetup->fwdSize(); 
         for(size_t i = 0; i < this->mDiff.size(); i++)
         {
            int cols = this->mspSetup->mult()(i);
            int specRows = this->mProj.at(i).cols();
            rSpecVal.block(0, start, specRows, cols) = this->mLl.bottomRows(specRows).asDiagonal()*this->mProj.at(i).transpose()*this->mWeights.asDiagonal()*physVal.block(0,start, physRows, cols);
            start += cols;
         }

      } else if(integrator == AssociatedLegendreTransform::IntegratorType::INTGDIVLL)
      { 
         // Compute integration
         int start = 0;
         int physRows = this->mspSetup->fwdSize(); 
         for(size_t i = 0; i < this->mDiff.size(); i++)
         {
            int cols = this->mspSetup->mult()(i);
            int specRows = this->mProj.at(i).cols();
            rSpecVal.block(0, start, specRows, cols) = this->mDivLl.bottomRows(specRows).asDiagonal()*this->mProj.at(i).transpose()*this->mWeights.asDiagonal()*physVal.block(0,start, physRows, cols);
            start += cols;
         }

      } else if(integrator == AssociatedLegendreTransform::IntegratorType::INTG)
      {
         // Compute integration
         int start = 0;
         int physRows = this->mspSetup->fwdSize(); 
         for(size_t i = 0; i < this->mDiff.size(); i++)
         {
            int cols = this->mspSetup->mult()(i);
            int specRows = this->mProj.at(i).cols();
            rSpecVal.block(0, start, specRows, cols) = this->mProj.at(i).transpose()*this->mWeights.asDiagonal()*physVal.block(0,start, physRows, cols);
            start += cols;
         }

      } else
      {
         throw Exception("Requested an unknown integrator");
      }
   }

   void AssociatedLegendreTransform::project(MatrixZ& rPhysVal, const MatrixZ& specVal, AssociatedLegendreTransform::ProjectorType::Id projector, Arithmetics::Id arithId)
   {
      assert(arithId == Arithmetics::SET);

      // assert right sizes for input  matrix
      assert(specVal.cols() == this->mspSetup->howmany());

      // assert right sizes for output matrix
      assert(rPhysVal.rows() == this->mspSetup->fwdSize());
      assert(rPhysVal.cols() == this->mspSetup->howmany());

      // Compute first derivative
      if(projector == AssociatedLegendreTransform::ProjectorType::DIFF)
      {
         int start = 0;
         int physRows = this->mspSetup->fwdSize(); 
         for(size_t i = 0; i < this->mDiff.size(); i++)
         {
            int cols = this->mspSetup->mult()(i);
            int specRows = this->mDiff.at(i).cols();
            rPhysVal.block(0, start, physRows, cols) = this->mDiff.at(i)*specVal.block(0,start, specRows, cols);
            start += cols;
         }

      // Compute \f$l(l+1)D\f$ projection
      } else if(projector == AssociatedLegendreTransform::ProjectorType::DIFFLL)
      {
         int start = 0;
         int physRows = this->mspSetup->fwdSize(); 
         for(size_t i = 0; i < this->mDivSin.size(); i++)
         {
            int cols = this->mspSetup->mult()(i);
            int specRows = this->mDivSin.at(i).cols();
            rPhysVal.block(0, start, physRows, cols) = this->mDiff.at(i)*this->mLl.bottomRows(specRows).asDiagonal()*specVal.block(0,start, specRows, cols);
            start += cols;
         }

      // Compute \f$l(l+1)/\sin\theta\f$ projection
      } else if(projector == AssociatedLegendreTransform::ProjectorType::DIVSIN)
      {
         int start = 0;
         int physRows = this->mspSetup->fwdSize(); 
         for(size_t i = 0; i < this->mDivSin.size(); i++)
         {
            int cols = this->mspSetup->mult()(i);
            int specRows = this->mDivSin.at(i).cols();
            rPhysVal.block(0, start, physRows, cols) = this->mDivSin.at(i)*specVal.block(0,start, specRows, cols);
            start += cols;
         }

      // Compute \f$1/\sin\theta\f$ projection
      } else if(projector == AssociatedLegendreTransform::ProjectorType::DIVSINLL)
      {
         int start = 0;
         int physRows = this->mspSetup->fwdSize(); 
         for(size_t i = 0; i < this->mDivSin.size(); i++)
         {
            int cols = this->mspSetup->mult()(i);
            int specRows = this->mDivSin.at(i).cols();
            rPhysVal.block(0, start, physRows, cols) = this->mDivSin.at(i)*this->mLl.bottomRows(specRows).asDiagonal()*specVal.block(0,start, specRows, cols);
            start += cols;
         }

      // Compute \f$1/\sin\theta \partial \sin\theta\f$ projection
      } else if(projector == AssociatedLegendreTransform::ProjectorType::DIVSINDIFFSIN)
      {
         throw Exception("DIVSINDIFFSIN not yet implemented");

      // Compute \f$l(l+1)\f$ projection
      } else if(projector == AssociatedLegendreTransform::ProjectorType::PROJLL)
      {
         int start = 0;
         int physRows = this->mspSetup->fwdSize(); 
         for(size_t i = 0; i < this->mProj.size(); i++)
         {
            int cols = this->mspSetup->mult()(i);
            int specRows = this->mProj.at(i).cols();
            rPhysVal.block(0, start, physRows, cols) = this->mProj.at(i)*this->mLl.bottomRows(specRows).asDiagonal()*specVal.block(0,start, specRows, cols);
            start += cols;
         }

      // Compute simple projection
      } else if(projector == AssociatedLegendreTransform::ProjectorType::PROJ)
      {
         int start = 0;
         int physRows = this->mspSetup->fwdSize(); 
         for(size_t i = 0; i < this->mProj.size(); i++)
         {
            int cols = this->mspSetup->mult()(i);
            int specRows = this->mProj.at(i).cols();
            rPhysVal.block(0, start, physRows, cols) = this->mProj.at(i)*specVal.block(0,start, specRows, cols);
            start += cols;
         }

      } else
      {
         throw Exception("Requested an unknown projector");
      }
   }

#ifdef GEOMHDISCC_STORAGEPROFILE
   MHDFloat AssociatedLegendreTransform::requiredStorage() const
   {
      MHDFloat mem = 0.0;

      // Storage for the grid and weight
      mem += static_cast<MHDFloat>(Debug::MemorySize<MHDFloat>::BYTES)*this->mXGrid.size();
      mem += static_cast<MHDFloat>(Debug::MemorySize<MHDFloat>::BYTES)*this->mThGrid.size();
      mem += static_cast<MHDFloat>(Debug::MemorySize<MHDFloat>::BYTES)*this->mWeights.size();

      // Storage for the projector
      for(size_t i = 0; i < this->mProj.size(); i++)
      {
         mem += static_cast<MHDFloat>(Debug::MemorySize<MHDFloat>::BYTES)*this->mProj.at(i).rows()*this->mProj.at(i).cols();
      }

      // Storage for the derivative
      for(size_t i = 0; i < this->mDiff.size(); i++)
      {
         mem += static_cast<MHDFloat>(Debug::MemorySize<MHDFloat>::BYTES)*this->mDiff.at(i).rows()*this->mDiff.at(i).cols();
      }

      // Storage for the 1/\sin\theta
      for(size_t i = 0; i < this->mDivSin.size(); i++)
      {
         mem += static_cast<MHDFloat>(Debug::MemorySize<MHDFloat>::BYTES)*this->mDivSin.at(i).rows()*this->mDivSin.at(i).cols();
      }

      return mem;
   }
#endif // GEOMHDISCC_STORAGEPROFILE

}
}
