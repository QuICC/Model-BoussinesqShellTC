/** 
 * @file DiagnosticCoordinator.cpp
 * @brief Source of the diagnostic coordinator
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
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
#if defined GEOMHDISCC_SPATIALSCHEME_FFF || defined GEOMHDISCC_SPATIALSCHEME_TFF || defined GEOMHDISCC_SPATIALSCHEME_TFT || defined GEOMHDISCC_SPATIALSCHEME_TTT 
   #include "Diagnostics/StreamVerticalWrapper.hpp"
   #include "Diagnostics/CartesianCflWrapper.hpp"
#elif defined GEOMHDISCC_SPATIALSCHEME_SLFL_TORPOL || defined GEOMHDISCC_SPATIALSCHEME_SLFM_TORPOL || defined GEOMHDISCC_SPATIALSCHEME_BLFL_TORPOL || defined GEOMHDISCC_SPATIALSCHEME_BLFMTORPOL
   #include "Diagnostics/SphericalTorPolWrapper.hpp"
   #include "Diagnostics/SphericalCflWrapper.hpp"
#endif //defined GEOMHDISCC_SPATIALSCHEME_FFF || defined GEOMHDISCC_SPATIALSCHEME_TFF || defined GEOMHDISCC_SPATIALSCHEME_TFT || defined GEOMHDISCC_SPATIALSCHEME_TTT 

namespace GeoMHDiSCC {

namespace Diagnostics {

   DiagnosticCoordinator::DiagnosticCoordinator()
      : mcMaxStep(0.1), mcMinStep(1e-11), mFixedStep(-1), mMaxError(-1.0), mCfl(0.0), mStartTime(0.0), mStartTimestep(0.0)
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
         this->mFixedStep = std::max(this->mcMinStep, tstep(1));
         this->mFixedStep = std::min(this->mFixedStep, this->mcMaxStep);

      #if defined GEOMHDISCC_SPATIALSCHEME_FFF || defined GEOMHDISCC_SPATIALSCHEME_TFF || defined GEOMHDISCC_SPATIALSCHEME_TFT || defined GEOMHDISCC_SPATIALSCHEME_TTT 
      // Create a stream function and vertical velocity wrapper
      } else if(scalars.count(PhysicalNames::STREAMFUNCTION) && scalars.count(PhysicalNames::VELOCITYZ))
      {
         SharedStreamVerticalWrapper spVelocity = SharedStreamVerticalWrapper(new StreamVerticalWrapper(scalars.find(PhysicalNames::STREAMFUNCTION)->second, scalars.find(PhysicalNames::VELOCITYZ)->second));

         this->mspCflWrapper = SharedCartesianCflWrapper(new CartesianCflWrapper(spVelocity, mesh));

         this->mFixedStep = tstep(1);

      #elif defined GEOMHDISCC_SPATIALSCHEME_SLFL_TORPOL || defined GEOMHDISCC_SPATIALSCHEME_SLFM_TORPOL || defined GEOMHDISCC_SPATIALSCHEME_BLFL_TORPOL || defined GEOMHDISCC_SPATIALSCHEME_BLFMTORPOL
      // Create a stream function and vertical velocity wrapper
      } else if(vectors.count(PhysicalNames::VELOCITY))
      {
         SharedSphericalTorPolWrapper spVelocity = SharedSphericalTorPolWrapper(new SphericalTorPolWrapper(vectors.find(PhysicalNames::VELOCITY)->second));

         this->mspCflWrapper = SharedSphericalCflWrapper(new SphericalCflWrapper(spVelocity, mesh));

         this->mFixedStep = tstep(1);
      #endif //defined GEOMHDISCC_SPATIALSCHEME_FFF || defined GEOMHDISCC_SPATIALSCHEME_TFF || defined GEOMHDISCC_SPATIALSCHEME_TFT || defined GEOMHDISCC_SPATIALSCHEME_TTT 

      // Required wrapper is not implemented
      } else
      {
         this->mFixedStep = tstep(1);
      }

      // Store configuration file start time
      this->mStartTime = tstep(0);

      // Store configuration file start step
      this->mStartTimestep = tstep(1);

      // Store error goal from configuration file (not enabled if fixed timestep is used)
      if(tstep(2) > 0)
      {
         this->mMaxError = tstep(2);
      }
   }

   void DiagnosticCoordinator::initialCfl()
   {
      // Used fixed timestep
      if(this->mFixedStep > 0 && this->mMaxError > 0)
      {
         this->mCfl = this->mcMinStep;

      } else if(this->mFixedStep > 0)
      {
         this->mCfl = this->mFixedStep;

      // Compute initial CFL condition
      } else if(this->mspCflWrapper)
      {
         // Compute CFL for initial state
         this->mCfl = this->mspCflWrapper->initialCfl();

         this->mCfl = std::min(this->mcMinStep, this->mCfl);
      }
   }

   void DiagnosticCoordinator::updateCfl()
   {
      // Used fixed timestep
      if(this->mFixedStep > 0)
      {
         this->mCfl = this->mFixedStep;

      // Compute CFL condition
      } else if(this->mspCflWrapper)
      {
         // Safety assert
         assert(this->mspCflWrapper);

         this->mCfl = this->mspCflWrapper->cfl();

         // Check for maximum timestep
         this->mCfl = std::min(this->mcMaxStep, this->mCfl);

         this->mCfl = std::min(-this->mFixedStep, this->mCfl);
      }
   }

   void DiagnosticCoordinator::synchronize()
   {
      //
      // Start of MPI block
      //
      #ifdef GEOMHDISCC_MPI

      if(this->mFixedStep <= 0 && this->mspCflWrapper)
      {
         // Reduce CFL on all CPUs to the global minimum
         MPI_Allreduce(MPI_IN_PLACE, &this->mCfl, 1, Parallel::MpiTypes::type<MHDFloat>(), MPI_MIN, MPI_COMM_WORLD);
      }

      //
      // End of MPI block
      //
      #endif // GEOMHDISCC_MPI
   }

   MHDFloat DiagnosticCoordinator::maxError() const
   {
      return this->mMaxError;
   }

   MHDFloat DiagnosticCoordinator::cfl() const
   {
      return this->mCfl;
   }

   MHDFloat DiagnosticCoordinator::startTime() const
   {
      return this->mStartTime;
   }

   MHDFloat DiagnosticCoordinator::startTimestep() const
   {
      return this->mStartTimestep;
   }

   void DiagnosticCoordinator::useStateTime(const MHDFloat time, const MHDFloat timestep)
   {
      // Configuration requests use of state time
      if(this->mStartTime < 0)
      {
         this->mStartTime = time;
      }

      // Configuration requests use of state timestep
      if(this->mStartTimestep < 0)
      {
         this->mStartTimestep = timestep;
      }
   }

}
}
