/** \file DiagnosticCoordinator.cpp
 *  \brief Source of the diagnostic coordinator
 */

// Debug includes
//
#include <cassert>
#include "Debug/DebuggerMacro.h"

// System includes
//

// External includes
//

// Class include
//
#include "Diagnostics/DiagnosticCoordinator.hpp"

// Project includes
//
#include "Base/MpiTypes.hpp"
#include "Exceptions/Exception.hpp"
#include "Diagnostics/StreamVerticalWrapper.hpp"

namespace GeoMHDiSCC {

namespace Diagnostics {

   DiagnosticCoordinator::DiagnosticCoordinator()
      : mcCourant(0.65), mFixedStep(-1), mCfl(0.0), mKinetic(0.0)
   {
   }

   DiagnosticCoordinator::~DiagnosticCoordinator()
   {
   }

   void DiagnosticCoordinator::init(const std::vector<Array>& mesh, const std::map<PhysicalNames::Id, Datatypes::SharedScalarVariableType>&  scalars, const std::map<PhysicalNames::Id, Datatypes::SharedVectorVariableType>&  vectors, const Array& tstep)
   {
      // Check for constant timestep setup
      if(tstep(1) > 0)
      {
         this->mFixedStep = tstep(1);

      // Create a stream function and vertical velocity wrapper
      } else if(scalars.count(PhysicalNames::STREAMFUNCTION) && scalars.count(PhysicalNames::VELOCITYZ))
      {
         SharedStreamVerticalWrapper   spWrapper = SharedStreamVerticalWrapper(new StreamVerticalWrapper(scalars.find(PhysicalNames::STREAMFUNCTION)->second, scalars.find(PhysicalNames::VELOCITYZ)->second));

         this->mspVelocityWrapper = spWrapper;

      // Required wrapper is not implemented
      } else
      {
         throw Exception("Unknown velocity wrapper is required!");
      }

      // Prepare mesh if CFL computation is required
      if(this->mFixedStep <= 0)
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
   }

   void DiagnosticCoordinator::initialCfl()
   {
      // Used fixed timestep
      if(this->mFixedStep > 0)
      {
         this->mCfl = this->mFixedStep;

      // Compute initial CFL condition
      } else
      {
         assert(this->mspVelocityWrapper);

         // Compute CFL for initial state
         this->updateCfl();

         // Assume a velocity of 100 to avoid problems with "zero" starting values 
         this->mCfl = std::min(this->mCfl, this->mcCourant*this->mMeshSpacings.at(0).minCoeff()/100.);
         this->mCfl = std::min(this->mCfl, this->mcCourant*this->mMeshSpacings.at(1).minCoeff()/100.);
         this->mCfl = std::min(this->mCfl, this->mcCourant*this->mMeshSpacings.at(2).minCoeff()/100.);
      }
   }

   void DiagnosticCoordinator::updateCfl()
   {
      // Used fixed timestep
      if(this->mFixedStep > 0)
      {
         this->mCfl = this->mFixedStep;

      // Compute initial CFL condition
      } else
      {
         // Safety assert
         assert(this->mspVelocityWrapper);

         // Compute most stringent CFL condition
         this->mCfl = this->mcCourant*this->mMeshSpacings.at(0).minCoeff()/this->mspVelocityWrapper->one().data().array().abs().maxCoeff();
         DebuggerMacro_showValue("Raw CFL Dx = ", 2, this->mMeshSpacings.at(0).minCoeff());
         DebuggerMacro_showValue("Raw CFL Vx = ", 2, this->mspVelocityWrapper->one().data().array().abs().maxCoeff());
         this->mCfl = std::min(this->mCfl, this->mcCourant*this->mMeshSpacings.at(1).minCoeff()/this->mspVelocityWrapper->two().data().array().abs().maxCoeff());
         DebuggerMacro_showValue("Raw CFL Dy = ", 2, this->mMeshSpacings.at(1).minCoeff());
         DebuggerMacro_showValue("Raw CFL Vy = ", 2, this->mspVelocityWrapper->two().data().array().abs().maxCoeff());
         this->mCfl = std::min(this->mCfl, this->mcCourant*this->mMeshSpacings.at(2).minCoeff()/this->mspVelocityWrapper->three().data().array().abs().maxCoeff());
         DebuggerMacro_showValue("Raw CFL Dz = ", 2, this->mMeshSpacings.at(2).minCoeff());
         DebuggerMacro_showValue("Raw CFL Vz = ", 2, this->mspVelocityWrapper->three().data().array().abs().maxCoeff());

         DebuggerMacro_showValue("Raw CFL cfl = ", 2, this->mCfl);
         /// Compute CFL condition : \f$\alpha\frac{\Delta x}{|v_{max}|}\f$
         /// \mhdBug This is not a clean CFL but min(0.1,cfl)
         //this->mCfl = std::min(this->mcCourant, this->mCfl);
         this->mCfl = std::min(0.1, this->mCfl);
      }
   }

   void DiagnosticCoordinator::updateKineticEnergy()
   {
   }

   void DiagnosticCoordinator::synchronize()
   {

      //
      // Start of MPI block
      //
      #ifdef GEOMHDISCC_MPI

      if(this->mFixedStep <= 0)
      {
         // Reduce CFL on all CPUs to the global minimum
         MPI_Allreduce(MPI_IN_PLACE, &this->mCfl, 1, Parallel::MpiTypes::type<MHDFloat>(), MPI_MIN, MPI_COMM_WORLD);
      }

      //
      // End of MPI block
      //
      #endif // GEOMHDISCC_MPI
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
