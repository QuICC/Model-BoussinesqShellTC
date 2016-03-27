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

      // Initialise the quadrature grid and weights and operators
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
      this->mProjOp.insert(std::make_pair(ProjectorType::PROJ,std::vector<Matrix>()));
      this->mProjOp.find(ProjectorType::PROJ)->second.reserve(this->mspSetup->slow().size());
      this->mProjOp.insert(std::make_pair(ProjectorType::DIFF,std::vector<Matrix>()));
      this->mProjOp.find(ProjectorType::DIFF)->second.reserve(this->mspSetup->slow().size());
      this->mProjOp.insert(std::make_pair(ProjectorType::DIVSIN,std::vector<Matrix>()));
      this->mProjOp.find(ProjectorType::DIVSIN)->second.reserve(this->mspSetup->slow().size());

      #ifdef GEOMHDISCC_MEMORYUSAGE_HIGH
      // Reserve storage for the weighted projectors, 1/sin projectors and derivative
      this->mIntgOp.insert(std::make_pair(IntegratorType::INTG,std::vector<Matrix>()));
      this->mIntgOp.find(IntegratorType::INTG)->second.reserve(this->mspSetup->slow().size());
      this->mIntgOp.insert(std::make_pair(IntegratorType::INTGDIFF,std::vector<Matrix>()));
      this->mIntgOp.find(IntegratorType::INTGDIFF)->second.reserve(this->mspSetup->slow().size());
      this->mIntgOp.insert(std::make_pair(IntegratorType::INTGDIVSIN,std::vector<Matrix>()));
      this->mIntgOp.find(IntegratorType::INTGDIVSIN)->second.reserve(this->mspSetup->slow().size());
      #endif //GEOMHDISCC_MEMORYUSAGE_HIGH

      // Loop over harmonic orders
      Matrix op;
      internal::Matrix  itmp, ipoly;
      for(int iM = 0; iM < this->mspSetup->slow().size(); iM++)
      {
         int m = this->mspSetup->slow()(iM);
         int maxL = this->mspSetup->fast().at(iM)(this->mspSetup->fast().at(iM).size()-1);

         // Allocate memory for the projector, 1/sin projector and derivatives
         this->mProjOp.find(ProjectorType::PROJ)->second.push_back(Matrix(this->mspSetup->fast().at(iM).size(), this->mThGrid.size()));
         this->mProjOp.find(ProjectorType::DIVSIN)->second.push_back(Matrix(this->mspSetup->fast().at(iM).size(), this->mThGrid.size()));
         this->mProjOp.find(ProjectorType::DIFF)->second.push_back(Matrix(this->mspSetup->fast().at(iM).size(), this->mThGrid.size()));

         op.resize(this->mThGrid.size(), maxL - m + 1);

         // Loop over harmonic degrees for projectors
         Polynomial::AssociatedLegendrePolynomial::Plm(op, ipoly, m, igrid);
         std::map<ProjectorType::Id,std::vector<Matrix> >::iterator projIt = this->mProjOp.find(ProjectorType::PROJ);
         for(int iL = 0; iL < this->mspSetup->fast().at(iM).size(); iL++)
         {
            int l = this->mspSetup->fast().at(iM)(iL);

            projIt->second.at(iM).row(iL) = op.col(l - m).transpose();
         }

         // Loop over harmonic degrees for derivative
         Polynomial::AssociatedLegendrePolynomial::dPlm(op, itmp, m, ipoly, igrid);
         projIt = this->mProjOp.find(ProjectorType::DIFF);
         for(int iL = 0; iL < this->mspSetup->fast().at(iM).size(); iL++)
         {
            int l = this->mspSetup->fast().at(iM)(iL);

            projIt->second.at(iM).row(iL) = op.col(l - m).transpose();
         }

         // Loop over harmonic degrees for 1/\sin\theta
         Polynomial::AssociatedLegendrePolynomial::sin_1Plm(op, itmp, m, igrid);
         projIt = this->mProjOp.find(ProjectorType::DIVSIN);
         for(int iL = 0; iL < this->mspSetup->fast().at(iM).size(); iL++)
         {
            int l = this->mspSetup->fast().at(iM)(iL);

            projIt->second.at(iM).row(iL) = op.col(l - m).transpose();
         }

         #ifdef GEOMHDISCC_MEMORYUSAGE_HIGH
         // Allocate memory for the weighted projector, 1/sin projector and derivatives
         this->mIntgOp.find(IntegratorType::INTG)->second.push_back(Matrix(this->mThGrid.size(), this->mspSetup->fast().at(iM).size()));
         this->mIntgOp.find(IntegratorType::INTGDIVSIN)->second.push_back(Matrix(this->mThGrid.size(), this->mspSetup->fast().at(iM).size()));
         this->mIntgOp.find(IntegratorType::INTGDIFF)->second.push_back(Matrix(this->mThGrid.size(), this->mspSetup->fast().at(iM).size()));

         this->mIntgOp.find(IntegratorType::INTG)->second.at(iM) = (this->mProjOp.find(ProjectorType::PROJ)->second.at(iM)*this->mWeights.asDiagonal()).transpose();
         this->mIntgOp.find(IntegratorType::INTGDIVSIN)->second.at(iM) = (this->mProjOp.find(ProjectorType::DIVSIN)->second.at(iM)*this->mWeights.asDiagonal()).transpose();
         this->mIntgOp.find(IntegratorType::INTGDIFF)->second.at(iM) = (this->mProjOp.find(ProjectorType::DIFF)->second.at(iM)*this->mWeights.asDiagonal()).transpose();
         #endif //GEOMHDISCC_MEMORYUSAGE_HIGH
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
         #ifdef GEOMHDISCC_MEMORYUSAGE_HIGH
         this->setIntegrator(rSpecVal, physVal, this->mIntgOp.find(integrator)->second);
         #else
         this->setIntegrator(rSpecVal, physVal, this->mProjOp.find(ProjectorType::DIFF)->second);
         #endif //GEOMHDISCC_MEMORYUSAGE_HIGH

      } else if(integrator == AssociatedLegendreTransform::IntegratorType::INTGLLDIFF)
      { 
         #ifdef GEOMHDISCC_MEMORYUSAGE_HIGH
         this->setMultLIntegrator(rSpecVal, physVal, this->mIntgOp.find(IntegratorType::INTGDIFF)->second, this->mLl);
         #else
         this->setMultLIntegrator(rSpecVal, physVal, this->mProjOp.find(ProjectorType::DIFF)->second, this->mLl);
         #endif //GEOMHDISCC_MEMORYUSAGE_HIGH

      } else if(integrator == AssociatedLegendreTransform::IntegratorType::INTGDIVLLDIFF)
      { 
         #ifdef GEOMHDISCC_MEMORYUSAGE_HIGH
         this->setMultLIntegrator(rSpecVal, physVal, this->mIntgOp.find(IntegratorType::INTGDIFF)->second, this->mDivLl);
         #else
         this->setMultLIntegrator(rSpecVal, physVal, this->mProjOp.find(ProjectorType::DIFF)->second, this->mDivLl);
         #endif //GEOMHDISCC_MEMORYUSAGE_HIGH

      } else if(integrator == AssociatedLegendreTransform::IntegratorType::INTGLLDIVSIN)
      { 
         #ifdef GEOMHDISCC_MEMORYUSAGE_HIGH
         this->setMultLIntegrator(rSpecVal, physVal, this->mIntgOp.find(IntegratorType::INTGDIVSIN)->second, this->mLl);
         #else
         this->setMultLIntegrator(rSpecVal, physVal, this->mProjOp.find(ProjectorType::DIVSIN)->second, this->mLl);
         #endif //GEOMHDISCC_MEMORYUSAGE_HIGH

      } else if(integrator == AssociatedLegendreTransform::IntegratorType::INTGDIVLLDIVSIN)
      { 
         #ifdef GEOMHDISCC_MEMORYUSAGE_HIGH
         this->setMultLIntegrator(rSpecVal, physVal, this->mIntgOp.find(IntegratorType::INTGDIVSIN)->second, this->mDivLl);
         #else
         this->setMultLIntegrator(rSpecVal, physVal, this->mProjOp.find(ProjectorType::DIVSIN)->second, this->mDivLl);
         #endif //GEOMHDISCC_MEMORYUSAGE_HIGH

      } else if(integrator == AssociatedLegendreTransform::IntegratorType::INTGDIVSIN)
      { 
         #ifdef GEOMHDISCC_MEMORYUSAGE_HIGH
         this->setIntegrator(rSpecVal, physVal, this->mIntgOp.find(integrator)->second);
         #else
         this->setIntegrator(rSpecVal, physVal, this->mProjOp.find(ProjectorType::DIVSIN)->second);
         #endif //GEOMHDISCC_MEMORYUSAGE_HIGH

      } else if(integrator == AssociatedLegendreTransform::IntegratorType::INTGLL2)
      { 
         #ifdef GEOMHDISCC_MEMORYUSAGE_HIGH
         this->setMultLIntegrator(rSpecVal, physVal, this->mIntgOp.find(IntegratorType::INTG)->second, this->mLl2);
         #else
         this->setMultLIntegrator(rSpecVal, physVal, this->mProjOp.find(ProjectorType::PROJ)->second, this->mLl2);
         #endif //GEOMHDISCC_MEMORYUSAGE_HIGH

      } else if(integrator == AssociatedLegendreTransform::IntegratorType::INTGLL)
      { 
         #ifdef GEOMHDISCC_MEMORYUSAGE_HIGH
         this->setMultLIntegrator(rSpecVal, physVal, this->mIntgOp.find(IntegratorType::INTG)->second, this->mLl);
         #else
         this->setMultLIntegrator(rSpecVal, physVal, this->mProjOp.find(ProjectorType::PROJ)->second, this->mLl);
         #endif //GEOMHDISCC_MEMORYUSAGE_HIGH

      } else if(integrator == AssociatedLegendreTransform::IntegratorType::INTGDIVLL)
      { 
         #ifdef GEOMHDISCC_MEMORYUSAGE_HIGH
         this->setMultLIntegrator(rSpecVal, physVal, this->mIntgOp.find(IntegratorType::INTG)->second, this->mDivLl);
         #else
         this->setMultLIntegrator(rSpecVal, physVal, this->mProjOp.find(ProjectorType::PROJ)->second, this->mDivLl);
         #endif //GEOMHDISCC_MEMORYUSAGE_HIGH

      } else if(integrator == AssociatedLegendreTransform::IntegratorType::INTG)
      {
         #ifdef GEOMHDISCC_MEMORYUSAGE_HIGH
         this->setIntegrator(rSpecVal, physVal, this->mIntgOp.find(integrator)->second);
         #else
         this->setIntegrator(rSpecVal, physVal, this->mProjOp.find(ProjectorType::PROJ)->second);
         #endif //GEOMHDISCC_MEMORYUSAGE_HIGH

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
         this->setProjector(rPhysVal, specVal, this->mProjOp.find(projector)->second);

      // Compute \f$l(l+1)D\f$ projection
      } else if(projector == AssociatedLegendreTransform::ProjectorType::DIFFLL)
      {
         this->setMultLProjector(rPhysVal, specVal, this->mProjOp.find(ProjectorType::DIFF)->second, this->mLl);

      // Compute \f$l(l+1)/\sin\theta\f$ projection
      } else if(projector == AssociatedLegendreTransform::ProjectorType::DIVSIN)
      {
         this->setProjector(rPhysVal, specVal, this->mProjOp.find(projector)->second);

      // Compute \f$1/\sin\theta\f$ projection
      } else if(projector == AssociatedLegendreTransform::ProjectorType::DIVSINLL)
      {
         this->setMultLProjector(rPhysVal, specVal, this->mProjOp.find(ProjectorType::DIVSIN)->second, this->mLl);

      // Compute \f$1/\sin\theta \partial \sin\theta\f$ projection
      } else if(projector == AssociatedLegendreTransform::ProjectorType::DIVSINDIFFSIN)
      {
         throw Exception("DIVSINDIFFSIN not yet implemented");

      // Compute \f$l(l+1)\f$ projection
      } else if(projector == AssociatedLegendreTransform::ProjectorType::PROJLL)
      {
         this->setMultLProjector(rPhysVal, specVal, this->mProjOp.find(ProjectorType::PROJ)->second, this->mLl);

      // Compute simple projection
      } else if(projector == AssociatedLegendreTransform::ProjectorType::PROJ)
      {
         this->setProjector(rPhysVal, specVal, this->mProjOp.find(projector)->second);

      } else
      {
         throw Exception("Requested an unknown projector");
      }
   }

   void AssociatedLegendreTransform::setIntegrator(MatrixZ& rSpecVal, const MatrixZ& physVal, const std::vector<Matrix>& ops)
   {
      // Compute integration
      int start = 0;
      int physRows = this->mspSetup->fwdSize(); 
      for(size_t i = 0; i < ops.size(); i++)
      {
         int cols = this->mspSetup->mult()(i);
         #ifdef GEOMHDISCC_MEMORYUSAGE_HIGH
         int specRows = ops.at(i).cols();
         rSpecVal.block(0, start, specRows, cols) = ops.at(i).transpose()*physVal.block(0,start, physRows, cols);
         #else
         int specRows = ops.at(i).rows();
         rSpecVal.block(0, start, specRows, cols) = ops.at(i)*(this->mWeights.asDiagonal()*physVal.block(0,start, physRows, cols));
         #endif //GEOMHDISCC_MEMORYUSAGE_HIGH
         start += cols;
      }
   }

   void AssociatedLegendreTransform::setMultLIntegrator(MatrixZ& rSpecVal, const MatrixZ& physVal, const std::vector<Matrix>& ops, const Array& mult)
   {
      // Compute integration
      int start = 0;
      int physRows = this->mspSetup->fwdSize(); 
      for(size_t i = 0; i < ops.size(); i++)
      {
         int cols = this->mspSetup->mult()(i);
         #ifdef GEOMHDISCC_MEMORYUSAGE_HIGH
         int specRows = ops.at(i).cols();
         rSpecVal.block(0, start, specRows, cols) = mult.bottomRows(specRows).asDiagonal()*(ops.at(i).transpose()*physVal.block(0,start, physRows, cols));
         #else
         int specRows = ops.at(i).rows();
         rSpecVal.block(0, start, specRows, cols) = mult.bottomRows(specRows).asDiagonal()*(ops.at(i)*(this->mWeights.asDiagonal()*physVal.block(0,start, physRows, cols)));
         #endif //GEOMHDISCC_MEMORYUSAGE_HIGH
         start += cols;
      }
   }

   void AssociatedLegendreTransform::setProjector(MatrixZ& rPhysVal, const MatrixZ& specVal, const std::vector<Matrix>& ops)
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

   void AssociatedLegendreTransform::setMultLProjector(MatrixZ& rPhysVal, const MatrixZ& specVal, const std::vector<Matrix>& ops, const Array& mult)
   {
      int start = 0;
      int physRows = this->mspSetup->fwdSize(); 
      for(size_t i = 0; i < ops.size(); i++)
      {
         int cols = this->mspSetup->mult()(i);
         int specRows = ops.at(i).rows();
         rPhysVal.block(0, start, physRows, cols) = ops.at(i).transpose()*(mult.bottomRows(specRows).asDiagonal()*specVal.block(0,start, specRows, cols));
         start += cols;
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
      std::map<ProjectorType::Id,std::vector<Matrix> >::const_iterator projIt = this->mProjOp.find(ProjectorType::PROJ);
      for(size_t i = 0; i < projIt->second.size(); i++)
      {
         mem += static_cast<MHDFloat>(Debug::MemorySize<MHDFloat>::BYTES)*projIt->second.at(i).rows()*projIt->second.at(i).cols();
      }

      // Storage for the derivative
      projIt = this->mProjOp.find(ProjectorType::DIFF);
      for(size_t i = 0; i < projIt->second.size(); i++)
      {
         mem += static_cast<MHDFloat>(Debug::MemorySize<MHDFloat>::BYTES)*projIt->second.at(i).rows()*projIt->second.at(i).cols();
      }

      // Storage for the 1/\sin\theta
      projIt = this->mProjOp.find(ProjectorType::DIVSIN);
      for(size_t i = 0; i < projIt->second.size(); i++)
      {
         mem += static_cast<MHDFloat>(Debug::MemorySize<MHDFloat>::BYTES)*projIt->second.at(i).rows()*projIt->second.at(i).cols();
      }

      #ifdef GEOMHDISCC_MEMORYUSAGE_HIGH
      // Storage for the integrator
      std::map<IntegratorType::Id,std::vector<Matrix> >::const_iterator intgIt = this->mIntgOp.find(IntegratorType::INTG);
      for(size_t i = 0; i < intgIt->second.size(); i++)
      {
         mem += static_cast<MHDFloat>(Debug::MemorySize<MHDFloat>::BYTES)*intgIt->second.at(i).rows()*intgIt->second.at(i).cols();
      }

      // Storage for the derivative
      intgIt = this->mIntgOp.find(IntegratorType::INTGDIFF);
      for(size_t i = 0; i < intgIt->second.size(); i++)
      {
         mem += static_cast<MHDFloat>(Debug::MemorySize<MHDFloat>::BYTES)*intgIt->second.at(i).rows()*intgIt->second.at(i).cols();
      }

      // Storage for the 1/\sin\theta
      intgIt = this->mIntgOp.find(IntegratorType::INTGDIVSIN);
      for(size_t i = 0; i < intgIt->second.size(); i++)
      {
         mem += static_cast<MHDFloat>(Debug::MemorySize<MHDFloat>::BYTES)*intgIt->second.at(i).rows()*intgIt->second.at(i).cols();
      }
      #endif //GEOMHDISCC_MEMORYUSAGE_HIGH

      return mem;
   }
#endif // GEOMHDISCC_STORAGEPROFILE

}
}
