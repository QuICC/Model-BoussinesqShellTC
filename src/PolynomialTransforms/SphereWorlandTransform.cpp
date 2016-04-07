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

namespace GeoMHDiSCC {

namespace Transform {

   Array SphereWorlandTransform::generateGrid(const int size)
   {
      // Initialise grid storage
      Array grid(size);

      // Create Chebyshev grid
      for(int k = 0; k < size; k++)
      {
         grid(k) = std::cos((Math::PI)*(static_cast<MHDFloat>(k)+0.5)/static_cast<MHDFloat>(size));

         grid(k) = std::sqrt((grid(k)+1.)/2.);
      }

      return grid;
   }

   SphereWorlandTransform::SphereWorlandTransform()
   {
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
      this->mGrid.resize(this->mspSetup->fwdSize());
      this->mWeights.resize(this->mspSetup->fwdSize());

      // Set the grid
      this->mGrid = SphereWorlandTransform::generateGrid(this->mspSetup->fwdSize());

      // Set the weights
      this->mWeights.setConstant(Math::PI/this->mspSetup->fwdSize());

      // Reserve storage for the projectors, 1st derivative
      this->mProjOp.insert(std::make_pair(ProjectorType::PROJ,std::vector<Matrix>()));
      this->mProjOp.find(ProjectorType::PROJ)->second.reserve(this->mspSetup->slow().size());
      this->mProjOp.insert(std::make_pair(ProjectorType::DIFF,std::vector<Matrix>()));
      this->mProjOp.find(ProjectorType::DIFF)->second.reserve(this->mspSetup->slow().size());

      #ifdef GEOMHDISCC_MEMORYUSAGE_HIGH
      // Reserve storage for the weighted projectors 
      this->mIntgOp.insert(std::make_pair(IntegratorType::INTG,std::vector<Matrix>()));
      this->mIntgOp.find(IntegratorType::INTG)->second.reserve(this->mspSetup->slow().size());
      #endif //GEOMHDISCC_MEMORYUSAGE_HIGH
   }

   void SphereWorlandTransform::integrate(MatrixZ& rSpecVal, const MatrixZ& physVal, SphereWorlandTransform::IntegratorType::Id integrator, Arithmetics::Id arithId)
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

   void SphereWorlandTransform::project(MatrixZ& rPhysVal, const MatrixZ& specVal, SphereWorlandTransform::ProjectorType::Id projector, Arithmetics::Id arithId)
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

      // Compute simple projection
      } else if(projector == AssociatedLegendreTransform::ProjectorType::PROJ)
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

#ifdef GEOMHDISCC_STORAGEPROFILE
   MHDFloat SphereWorlandTransform::requiredStorage() const
   {
      MHDFloat mem = 0.0;

      // Storage for the grid and weight
      mem += static_cast<MHDFloat>(Debug::MemorySize<MHDFloat>::BYTES)*this->mGrid.size();
      mem += static_cast<MHDFloat>(Debug::MemorySize<MHDFloat>::BYTES)*this->mWeights.size();

      // Storage for the projector
      std::map<ProjectorType::Id,std::vector<Matrix> >::const_iterator projIt = this->mProjOp.find(ProjectorType::PROJ);
      for(size_t i = 0; i < projIt->second.size(); i++)
      {
         mem += static_cast<MHDFloat>(Debug::MemorySize<MHDFloat>::BYTES)*projIt->second.at(i).size();
      }

      // Storage for the derivative
      projIt = this->mProjOp.find(ProjectorType::DIFF);
      for(size_t i = 0; i < projIt->second.size(); i++)
      {
         mem += static_cast<MHDFloat>(Debug::MemorySize<MHDFloat>::BYTES)*projIt->second.at(i).size();
      }

      #ifdef GEOMHDISCC_MEMORYUSAGE_HIGH
      // Storage for the integrator
      std::map<IntegratorType::Id,std::vector<Matrix> >::const_iterator intgIt = this->mIntgOp.find(IntegratorType::INTG);
      for(size_t i = 0; i < intgIt->second.size(); i++)
      {
         mem += static_cast<MHDFloat>(Debug::MemorySize<MHDFloat>::BYTES)*intgIt->second.at(i).size();
      }
      #endif //GEOMHDISCC_MEMORYUSAGE_HIGH

      return mem;
   }
#endif // GEOMHDISCC_STORAGEPROFILE

}
}
