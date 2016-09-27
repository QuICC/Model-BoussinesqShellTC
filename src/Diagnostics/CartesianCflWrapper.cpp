/** 
 * @file CartesianCflWrapper.cpp
 * @brief Source of the CFL constraint wrapper in a Cartesian geometry
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Debug includes
//

// System includes
//

// External includes
//

// Class include
//
#include "Diagnostics/CartesianCflWrapper.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Diagnostics {

   CartesianCflWrapper::CartesianCflWrapper(const SharedIVelocityWrapper spVelocity, const std::vector<Array>& mesh)
      : ICflWrapper(spVelocity), mcCourant(0.65)
   {
      // Initialize the mesh
      this->initMesh(mesh);
   }

   CartesianCflWrapper::~CartesianCflWrapper()
   {
   }

   void CartesianCflWrapper::initMesh(const std::vector<Array>& mesh)
   {
      // Compute the mesh spacings
      this->mMeshSpacings.reserve(mesh.size());
      // Loop over all dimensions
      for(size_t i = 0; i < mesh.size(); i++)
      {
         // Create storage
         this->mMeshSpacings.push_back(Array(mesh.at(i).size()));

         // Extract minimal spacing for each grid in current direction
         for(int j = 0; j < mesh.at(i).size(); ++j)
         {
            // Get internal points grid spacing
            if(j > 0 && j < mesh.at(i).size() - 1)
            {
               this->mMeshSpacings.back()(j) = std::min(std::abs(mesh.at(i)(j) - mesh.at(i)(j-1)), std::abs(mesh.at(i)(j) - mesh.at(i)(j+1)));

            // Get left endpoint grid spacing
            } else if(j > 0)
            {
               this->mMeshSpacings.back()(j) = std::abs(mesh.at(i)(j) - mesh.at(i)(j-1));

            // Get right endpoint grid spacing
            } else
            {
               this->mMeshSpacings.back()(j) = std::abs(mesh.at(i)(j) - mesh.at(i)(j+1));
            }
         }
      }
   }

   MHDFloat CartesianCflWrapper::initialCfl() const
   {
      MHDFloat cfl = this->cfl();

      // Assume a velocity of 100 to avoid problems with "zero" starting values 
      cfl = std::min(cfl, this->mcCourant*this->mMeshSpacings.at(0).minCoeff()/100.);
      cfl = std::min(cfl, this->mcCourant*this->mMeshSpacings.at(1).minCoeff()/100.);
      cfl = std::min(cfl, this->mcCourant*this->mMeshSpacings.at(2).minCoeff()/100.);

      return cfl;
   }

   MHDFloat CartesianCflWrapper::cfl() const
   {
      // Compute most stringent CFL condition
      MHDFloat cfl = 1.0;

      int nK = this->mspVelocity->spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();

      // CFL from first component
      for(int k = 0; k < nK; ++k)
      {
         int k_ = this->mspVelocity->spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(k);
         cfl = std::min(cfl, this->mcCourant*this->mMeshSpacings.at(0)(k_)/this->mspVelocity->one().slice(k).array().abs().maxCoeff());
      }

      // CFL from second component
      for(int k = 0; k < nK; ++k)
      {
         int nJ = this->mspVelocity->spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(k);
         for(int j = 0; j < nJ; ++j)
         {
            int j_ = this->mspVelocity->spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT2D>(j,k);
            cfl = std::min(cfl, this->mcCourant*this->mMeshSpacings.at(1)(j_)/this->mspVelocity->two().profile(j,k).array().abs().maxCoeff());
         }
      }

      // CFL from second component
      for(int k = 0; k < nK; ++k)
      {
         int nJ = this->mspVelocity->spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(k);
         for(int j = 0; j < nJ; ++j)
         {
            int nI = this->mspVelocity->spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DATF1D>(k);
            for(int i = 0; i < nI; ++i)
            {
               cfl = std::min(cfl, this->mcCourant*this->mMeshSpacings.at(2)(i)/std::abs(this->mspVelocity->three().point(i,j,k)));
            }
         }
      }

      return cfl;
   }

}
}
