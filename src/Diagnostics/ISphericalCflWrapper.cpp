/** 
 * @file ISphericalCflWrapper.cpp
 * @brief Source of the CFL constraint wrapper in a spherical geometry
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
#include "Diagnostics/ISphericalCflWrapper.hpp"

// Project includes
//

namespace QuICC {

namespace Diagnostics {

   ISphericalCflWrapper::ISphericalCflWrapper(const SharedIVelocityWrapper spVelocity)
      : ICflWrapper(spVelocity), mcCourant(0.65)
   {
   }

   ISphericalCflWrapper::~ISphericalCflWrapper()
   {
   }

   void ISphericalCflWrapper::init(const std::vector<Array>& mesh)
   {
      // Initialize the mesh
      this->initMesh(mesh);
   }

   void ISphericalCflWrapper::initMesh(const std::vector<Array>& mesh)
   {
      // Compute the mesh spacings
      this->mMeshSpacings.reserve(2);

      // Storage for radial grid spacing
      this->mMeshSpacings.push_back(Array(mesh.at(0).size()));

      // Storage for horizontal average grid spacing
      this->mMeshSpacings.push_back(Array(mesh.at(0).size()));

      // Compute grid spacings
      for(int j = 0; j < mesh.at(0).size(); ++j)
      {
         // Get internal points grid spacing
         if(j > 0 && j < mesh.at(0).size() - 1)
         {
            this->mMeshSpacings.at(0)(j) = std::min(std::abs(mesh.at(0)(j) - mesh.at(0)(j-1)), std::abs(mesh.at(0)(j) - mesh.at(0)(j+1)));

         // Get left endpoint grid spacing
         } else if(j > 0)
         {
            this->mMeshSpacings.at(0)(j) = std::abs(mesh.at(0)(j) - mesh.at(0)(j-1));

            // Get right endpoint grid spacing
         } else
         {
            this->mMeshSpacings.at(0)(j) = std::abs(mesh.at(0)(j) - mesh.at(0)(j+1));
         }

         // Compute average horizontal grid spacing
         MHDFloat effL = this->effectiveMaxL(mesh.at(0)(j));
         this->mMeshSpacings.at(1)(j) = mesh.at(0)(j)/std::sqrt(effL*(effL + 1.0));
      }
   }

   MHDFloat ISphericalCflWrapper::initialCfl() const
   {
      MHDFloat cfl = this->cfl();

      // Assume a velocity of 100 to avoid problems with "zero" starting values 
      cfl = std::min(cfl, this->mcCourant*this->mMeshSpacings.at(0).minCoeff()/100.);
      cfl = std::min(cfl, this->mcCourant*this->mMeshSpacings.at(1).minCoeff()/100.);

      return cfl;
   }

   MHDFloat ISphericalCflWrapper::cfl() const
   {
      // Compute most stringent CFL condition
      MHDFloat cfl = 1.0;

      int nR = this->mspVelocity->spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();
      for(int i = 0; i < nR; ++i)
      {
         int iR = this->mspVelocity->spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(i);

         // Radial CFL
         cfl = std::min(cfl, this->mcCourant*this->mMeshSpacings.at(0)(iR)/this->mspVelocity->one().slice(i).array().abs().maxCoeff());
   
         // Horizontal CFL
         cfl = std::min(cfl, this->mcCourant*this->mMeshSpacings.at(1)(iR)/(this->mspVelocity->two().slice(i).array().pow(2) + this->mspVelocity->three().slice(i).array().pow(2)).array().sqrt().maxCoeff());
      }

      return cfl;
   }

}
}
