/**
 * @file SphereConserveAngularMomentumWriter.cpp
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
#include "IoVariable/SphereConserveAngularMomentumWriter.hpp"

// Project includes
//
#include "Enums/Dimensions.hpp"
#include "Enums/FieldIds.hpp"
#include "IoVariable/AngularMomentumTags.hpp"

namespace QuICC {

namespace IoVariable {

   SphereConserveAngularMomentumWriter::SphereConserveAngularMomentumWriter(const std::string& prefix, const std::string& type)
      : SphereAngularMomentumWriter(prefix, type), mUnitMomentum(3)
   {
   }

   SphereConserveAngularMomentumWriter::~SphereConserveAngularMomentumWriter()
   {
   }

   void SphereConserveAngularMomentumWriter::init()
   {
      this->mUnitMomentum(0) = -std::sqrt(128.0/75.0);
      this->mUnitMomentum(1) = std::sqrt(128.0/75.0);
      this->mUnitMomentum(2) = std::sqrt(64.0/75.0);

      SphereAngularMomentumWriter::init();
   }

   void SphereConserveAngularMomentumWriter::compute(Transform::TransformCoordinatorType& coord)
   {
      SphereAngularMomentumWriter::compute(coord);

      // get iterator to field
      vector_iterator vIt;
      vector_iterator_range vRange = this->vectorRange();
      assert(std::distance(vRange.first, vRange.second) == 1);
      assert(FieldComponents::Spectral::ONE == FieldComponents::Spectral::TOR);

      if(this->mHasMOrdering)
      {
         // Loop over harmonic order m
         for(int k = 0; k < this->res().cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
         {
            int m_ = this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k);

            for(int j = 0; j < this->res().cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); j++)
            {
               int l_ = this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j, k);

               if(l_== 1 && m_ == 0)
               {
                  MHDComplex p = vRange.first->second->rDom(0).rPerturbation().rComp(FieldComponents::Spectral::TOR).point(0,j,k);
                  vRange.first->second->rDom(0).rPerturbation().rComp(FieldComponents::Spectral::TOR).setPoint(p - this->mMomentum(2)/this->mUnitMomentum(2), 0,j,k);
               } else if(l_== 1 && m_ == 1)
               {
                  MHDComplex p = vRange.first->second->rDom(0).rPerturbation().rComp(FieldComponents::Spectral::TOR).point(0,j,k);
                  vRange.first->second->rDom(0).rPerturbation().rComp(FieldComponents::Spectral::TOR).setPoint(p - MHDComplex(this->mMomentum(0)/this->mUnitMomentum(0),this->mMomentum(1)/this->mUnitMomentum(1)), 0, j, k);
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
                  MHDComplex p = vRange.first->second->rDom(0).rPerturbation().rComp(FieldComponents::Spectral::TOR).point(0,j,k);
                  vRange.first->second->rDom(0).rPerturbation().rComp(FieldComponents::Spectral::TOR).setPoint(p - this->mMomentum(2)/this->mUnitMomentum(2), 0,j,k);
               } else if(l_ == 1 && m_ == 1)
               {
                  MHDComplex p = vRange.first->second->rDom(0).rPerturbation().rComp(FieldComponents::Spectral::TOR).point(0,j,k);
                  vRange.first->second->rDom(0).rPerturbation().rComp(FieldComponents::Spectral::TOR).setPoint(p - MHDComplex(this->mMomentum(0)/this->mUnitMomentum(0),this->mMomentum(1)/this->mUnitMomentum(1)), 0, j, k);
               }
            }
         }
      }
   }

}
}
