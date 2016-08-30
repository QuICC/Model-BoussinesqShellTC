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
      this->mGrid.resize(this->mspSetup->fwdSize());
      this->mWeights.resize(this->mspSetup->fwdSize());

      // Set the grid and weights
      internal::Array igrid, iweights;
      WorlandChebyshevRule::computeQuadrature(this->mGrid, this->mWeights, igrid, iweights, this->mspSetup->fwdSize());

      // Reserve storage for the projectors
      this->mProjOp.insert(std::make_pair(ProjectorType::PROJ,std::vector<Matrix>()));
      this->mProjOp.find(ProjectorType::PROJ)->second.reserve(this->mspSetup->slow().size());
      this->mProjOp.insert(std::make_pair(ProjectorType::DIVR,std::vector<Matrix>()));
      this->mProjOp.find(ProjectorType::DIVR)->second.reserve(this->mspSetup->slow().size());
      this->mProjOp.insert(std::make_pair(ProjectorType::DIFF,std::vector<Matrix>()));
      this->mProjOp.find(ProjectorType::DIFF)->second.reserve(this->mspSetup->slow().size());
      this->mProjOp.insert(std::make_pair(ProjectorType::DIVRDIFFR,std::vector<Matrix>()));
      this->mProjOp.find(ProjectorType::DIVRDIFFR)->second.reserve(this->mspSetup->slow().size());
      this->mProjOp.insert(std::make_pair(ProjectorType::LAPLH,std::vector<Matrix>()));
      this->mProjOp.find(ProjectorType::LAPLH)->second.reserve(this->mspSetup->slow().size());
      this->mProjOp.insert(std::make_pair(ProjectorType::DIVRLAPLH,std::vector<Matrix>()));
      this->mProjOp.find(ProjectorType::DIVRLAPLH)->second.reserve(this->mspSetup->slow().size());
      this->mProjOp.insert(std::make_pair(ProjectorType::DIFFLAPLH,std::vector<Matrix>()));
      this->mProjOp.find(ProjectorType::DIFFLAPLH)->second.reserve(this->mspSetup->slow().size());

      // Reserve storage for the weighted projectors 
      this->mIntgOp.insert(std::make_pair(IntegratorType::INTG,std::vector<Matrix>()));
      this->mIntgOp.find(IntegratorType::INTG)->second.reserve(this->mspSetup->slow().size());
      this->mIntgOp.insert(std::make_pair(IntegratorType::INTGI4DIVR,std::vector<Matrix>()));
      this->mIntgOp.find(IntegratorType::INTGI4DIVR)->second.reserve(this->mspSetup->slow().size());
      this->mIntgOp.insert(std::make_pair(IntegratorType::INTGI4DIVRDIFFR,std::vector<Matrix>()));
      this->mIntgOp.find(IntegratorType::INTGI4DIVRDIFFR)->second.reserve(this->mspSetup->slow().size());
      this->mIntgOp.insert(std::make_pair(IntegratorType::INTGI6DIVR,std::vector<Matrix>()));
      this->mIntgOp.find(IntegratorType::INTGI6DIVR)->second.reserve(this->mspSetup->slow().size());
      this->mIntgOp.insert(std::make_pair(IntegratorType::INTGI6DIVRDIFFR,std::vector<Matrix>()));
      this->mIntgOp.find(IntegratorType::INTGI6DIVRDIFFR)->second.reserve(this->mspSetup->slow().size());
      this->mIntgOp.insert(std::make_pair(IntegratorType::INTGI6LAPLH,std::vector<Matrix>()));
      this->mIntgOp.find(IntegratorType::INTGI6LAPLH)->second.reserve(this->mspSetup->slow().size());
      this->mIntgOp.insert(std::make_pair(IntegratorType::INTGINVLAPLH,std::vector<Matrix>()));
      this->mIntgOp.find(IntegratorType::INTGINVLAPLH)->second.reserve(this->mspSetup->slow().size());

      // Loop over harmonic degrees
      Matrix op;
      internal::Matrix  itmp, ipoly;
      for(int iL = 0; iL < this->mspSetup->slow().size(); iL++)
      {
         int l = this->mspSetup->slow()(iL);

         // Allocate memory for the projectors
         this->mProjOp.find(ProjectorType::PROJ)->second.push_back(Matrix(this->mspSetup->fast().at(iL).size(), this->mGrid.size()));
         this->mProjOp.find(ProjectorType::DIVR)->second.push_back(Matrix(this->mspSetup->fast().at(iL).size(), this->mGrid.size()));
         this->mProjOp.find(ProjectorType::DIFF)->second.push_back(Matrix(this->mspSetup->fast().at(iL).size(), this->mGrid.size()));
         this->mProjOp.find(ProjectorType::DIVRDIFFR)->second.push_back(Matrix(this->mspSetup->fast().at(iL).size(), this->mGrid.size()));
         this->mProjOp.find(ProjectorType::LAPLH)->second.push_back(Matrix(this->mspSetup->fast().at(iL).size(), this->mGrid.size()));
         this->mProjOp.find(ProjectorType::DIVRLAPLH)->second.push_back(Matrix(this->mspSetup->fast().at(iL).size(), this->mGrid.size()));
         this->mProjOp.find(ProjectorType::DIFFLAPLH)->second.push_back(Matrix(this->mspSetup->fast().at(iL).size(), this->mGrid.size()));

         op.resize(this->mGrid.size(), this->mspSetup->fast().at(iL).size());

         // Projector: P
         Polynomial::WorlandPolynomial::Wnl(op, ipoly, l, igrid);
         std::map<ProjectorType::Id,std::vector<Matrix> >::iterator projIt = this->mProjOp.find(ProjectorType::PROJ);
         projIt->second.at(iL) = op.transpose();

         // Allocate memory for the weighted integrator
         this->mIntgOp.find(IntegratorType::INTG)->second.push_back(Matrix(this->mGrid.size(), this->mspSetup->fast().at(iL).size()));
         this->mIntgOp.find(IntegratorType::INTG)->second.at(iL) = (this->mProjOp.find(ProjectorType::PROJ)->second.at(iL)*this->mWeights.asDiagonal()).transpose();
      }
   }

   void CylinderWorlandTransform::integrate(MatrixZ& rSpecVal, const MatrixZ& physVal, CylinderWorlandTransform::IntegratorType::Id integrator)
   {
      // assert right sizes for input matrix
      assert(physVal.rows() == this->mspSetup->fwdSize());
      assert(physVal.cols() == this->mspSetup->howmany());

      // assert right sizes for output matrix
      assert(rSpecVal.cols() == this->mspSetup->howmany());

      // Compute basic integrator
      if(integrator == CylinderWorlandTransform::IntegratorType::INTG)
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
      if(projector == CylinderWorlandTransform::ProjectorType::PROJ)
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

#ifdef GEOMHDISCC_STORAGEPROFILE
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
#endif // GEOMHDISCC_STORAGEPROFILE

}
}
