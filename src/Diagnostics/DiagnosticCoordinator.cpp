/** \file DiagnosticCoordinator.cpp
 *  \brief Source of the diagnostic coordinator
 */

// Debug includes
//
#include <cassert>

// System includes
//

// External includes
//

// Class include
//
#include "Diagnostics/DiagnosticCoordinator.hpp"

// Project includes
//
#include "Exceptions/Exception.hpp"
#include "Diagnostics/StreamVerticalWrapper.hpp"

#include <iostream>

namespace GeoMHDiSCC {

namespace Diagnostics {

   DiagnosticCoordinator::DiagnosticCoordinator()
      : mcCourant(0.65), mCfl(0.0), mKinetic(0.0)
   {
   }

   DiagnosticCoordinator::~DiagnosticCoordinator()
   {
   }

   void DiagnosticCoordinator::init(const std::vector<Array>& mesh, const std::map<PhysicalNames::Id, Datatypes::SharedScalarVariableType>&  scalars, const std::map<PhysicalNames::Id, Datatypes::SharedVectorVariableType>&  vectors)
   {
      // Create a stream function and vertical velocity wrapper
      if(scalars.count(PhysicalNames::STREAMFUNCTION) && scalars.count(PhysicalNames::VELOCITYZ))
      {
         SharedStreamVerticalWrapper   spWrapper = SharedStreamVerticalWrapper(new StreamVerticalWrapper(scalars.find(PhysicalNames::STREAMFUNCTION)->second, scalars.find(PhysicalNames::VELOCITYZ)->second));

         this->mspVelocityWrapper = spWrapper;

      // Required wrapper is not implemented
      } else
      {
         throw Exception("Unknown velocity wrapper is required!");
      }

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

   void DiagnosticCoordinator::updateCfl()
   {
      // Safety assert
      assert(this->mspVelocityWrapper);

      // Compute most stringent CFL condition
      this->mCfl = this->mcCourant*this->mMeshSpacings.at(0).minCoeff()/this->mspVelocityWrapper->one().data().array().abs().maxCoeff();
      this->mCfl = std::min(this->mCfl, this->mcCourant*this->mMeshSpacings.at(1).minCoeff()/this->mspVelocityWrapper->two().data().array().abs().maxCoeff());
      this->mCfl = std::min(this->mCfl, this->mcCourant*this->mMeshSpacings.at(2).minCoeff()/this->mspVelocityWrapper->three().data().array().abs().maxCoeff());

      std::cerr << "Raw CFL: " << this->mCfl << std::endl;
      /// Compute CFL condition : \f$\alpha\frac{\Delta x}{|v_{max}|}\f$
      //this->mCfl = std::min(this->mcCourant, this->mCfl);
      this->mCfl = std::min(0.1, this->mCfl);
   }

   void DiagnosticCoordinator::updateKineticEnergy()
   {
      // Safety assert
      assert(this->mspVelocityWrapper);

   }

   MHDFloat DiagnosticCoordinator::cfl() const
   {
      return this->mCfl;
   }

   MHDFloat DiagnosticCoordinator::kineticEnergy() const
   {
      return this->mKinetic;
   }

}
}
