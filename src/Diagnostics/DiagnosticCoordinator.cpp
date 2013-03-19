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

namespace GeoMHDiSCC {

namespace Diagnostics {

   DiagnosticCoordinator::DiagnosticCoordinator()
      : mcCourant(0.65), mCfl(0.0), mKinetic(0.0), mMinSpacing(0.0)
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

      // Get the minimal grid spacing
      /// \mhdBug Computationg of the minimal grid spacing is very strict
      Array minSpace(mesh.size());
      for(size_t i = 0; i < mesh.size(); i++)
      {
         minSpace(i) = (mesh.at(i).segment(0, mesh.at(i).size()-1) - mesh.at(i).segment(1, mesh.at(i).size()-1)).array().abs().minCoeff();
      }

      // Set the minimal spacing
      this->mMinSpacing = minSpace.array().minCoeff();
   }

   void DiagnosticCoordinator::updateCfl()
   {
      // Safety assert
      assert(this->mspVelocityWrapper);

      // Compute the maximum velocity amplitude on the grid
      MHDFloat maxVel = std::sqrt((this->mspVelocityWrapper->one().data().array().pow(2) + this->mspVelocityWrapper->two().data().array().pow(2) + this->mspVelocityWrapper->three().data().array().pow(2)).maxCoeff());

      /// Compute CFL condition : \f$\alpha\frac{\Delta x}{|v_{max}|}\f$
      this->mCfl = std::min(this->mcCourant, this->mcCourant*this->mMinSpacing/maxVel);
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
