/**
 * @file SphereAngularMomentumWriter.cpp
 * @brief Source of the implementation of the ASCII sphere angular momentum
 */

// Configuration includes
//
#include "Framework/FrameworkMacro.h"

// System includes
//
#include <iomanip>

// External includes
//

// Class include
//
#include "IoVariable/SphereAngularMomentumWriter.hpp"

// Project includes
//
#include "Enums/Dimensions.hpp"
#include "Enums/FieldIds.hpp"
#include "IoVariable/AngularMomentumTags.hpp"
#include "PolynomialTransforms/WorlandOperators.hpp"

namespace QuICC {

namespace IoVariable {

   SphereAngularMomentumWriter::SphereAngularMomentumWriter(const std::string& prefix, const std::string& type)
      : IVariableAsciiWriter(prefix + AngularMomentumTags::BASENAME, AngularMomentumTags::EXTENSION, prefix + AngularMomentumTags::HEADER, type, AngularMomentumTags::VERSION, Dimensions::Space::SPECTRAL, EXTEND), mMomentum(3)
   {
   }

   SphereAngularMomentumWriter::~SphereAngularMomentumWriter()
   {
   }

   void SphereAngularMomentumWriter::init()
   {
      #if defined QUICC_SPATIALSCHEME_BLFM || defined QUICC_SPATIALSCHEME_WLFM
         this->mHasMOrdering = true;
      #endif //defined QUICC_SPATIALSCHEME_BLFM || defined QUICC_SPATIALSCHEME_WLFM
      #if defined QUICC_SPATIALSCHEME_BLFL || defined QUICC_SPATIALSCHEME_WLFL
         this->mHasMOrdering = false;
      #endif // defined QUICC_SPATIALSCHEME_BLFL || defined QUICC_SPATIALSCHEME_WLFL

      IVariableAsciiWriter::init();

      this->mAngMomLM = MatrixI::Constant(2, 2, -1);

      if(this->mHasMOrdering)
      {
         // Loop over harmonic order m
         for(int k = 0; k < this->res().cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
         {
            int m_ = this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k);
            for(int j = 0; j < this->res().cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); j++)
            {
               int l_ = this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j, k);

               if(l_ == 1 && m_ == 0)
               {
                  this->mAngMomLM(0,0) = k;
                  this->mAngMomLM(1,0) = j;
               } else if(l_ == 1 && m_ == 1)
               {
                  this->mAngMomLM(0,1) = k;
                  this->mAngMomLM(1,1) = j;
               }
            }
         }
      } else
      {
         // Loop over harmonic degree l
         for(int k = 0; k < this->res().cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
         {
            int l_ = this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k);
            for(int j = 0; j < this->res().cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); j++)
            {
               int m_ = this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j,k);
               if(l_ == 1 && m_ == 0)
               {
                  this->mAngMomLM(0,0) = k;
                  this->mAngMomLM(1,0) = j;
               } else if(l_ == 1 && m_ == 1)
               {
                  this->mAngMomLM(0,1) = k;
                  this->mAngMomLM(1,1) = j;
               }
            }
         }
      }

      // Compute operator if required
      if(this->mAngMomLM.sum() > -4)
      {
         int nN = this->res().sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
         Polynomial::WorlandOperators::integralR3(this->mAngMomOp, nN, 1);
         assert(this->mAngMomOp.rows() == nN && this->mAngMomOp.cols() == 1);
      }
   }

   void SphereAngularMomentumWriter::compute(Transform::TransformCoordinatorType& coord)
   {
      // get iterator to field
      vector_iterator_range vRange = this->vectorRange();
      assert(std::distance(vRange.first, vRange.second) == 1);
      assert(FieldComponents::Spectral::ONE == FieldComponents::Spectral::TOR);

      ArrayZ mom;
      this->mMomentum.setZero();
      if(this->mAngMomLM(0,1) != -1)
      {
         MHDFloat c = std::sqrt(32.0*Math::PI/3.0);
         mom = (this->mAngMomOp.transpose()*vRange.first->second->dom(0).total().comp(FieldComponents::Spectral::TOR).profile(this->mAngMomLM(1,1),this->mAngMomLM(0,1)));
         this->mMomentum(0) = -c*mom(0).real();
         this->mMomentum(1) = c*mom(0).imag();
      }

      if(this->mAngMomLM(0,0) != -1)
      {
         MHDFloat c = std::sqrt(16.0*Math::PI/3.0);
         mom = (this->mAngMomOp.transpose()*vRange.first->second->dom(0).total().comp(FieldComponents::Spectral::TOR).profile(this->mAngMomLM(1,0),this->mAngMomLM(0,0)));
         this->mMomentum(2) = c*mom(0).real();
      }
   }

   void SphereAngularMomentumWriter::write()
   {
      // Create file
      this->preWrite();

      // Get the "global" value from MPI code
      #ifdef QUICC_MPI
         MPI_Allreduce(MPI_IN_PLACE, this->mMomentum.data(), this->mMomentum.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      #endif //QUICC_MPI

      // Check if the workflow allows IO to be performed
      if(FrameworkMacro::allowsIO())
      {
         this->mFile << std::setprecision(14) << this->mTime << "\t" << this->mMomentum.norm() << "\t" << this->mMomentum.transpose();
         
         // End line
         this->mFile << std::endl;
      }

      // Close file
      this->postWrite();

      // Abort if is NaN
      if(std::isnan(this->mMomentum.norm()))
      {
         FrameworkMacro::abort(99);

         throw std::logic_error("Sphere angular momentum is NaN!");
      }
   }

}
}
