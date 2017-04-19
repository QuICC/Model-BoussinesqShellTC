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
#if defined QUICC_SPATIALSCHEME_TFF_TORPOL
   #include "Diagnostics/CartesianTorPolWrapper.hpp"
   #include "Diagnostics/CartesianCflWrapper.hpp"
#elif defined QUICC_SPATIALSCHEME_SLFL_TORPOL || defined QUICC_SPATIALSCHEME_SLFM_TORPOL
   #include "Diagnostics/SphericalTorPolWrapper.hpp"
   #include "Diagnostics/ShellCflWrapper.hpp"
#elif defined QUICC_SPATIALSCHEME_BLFL_TORPOL || defined QUICC_SPATIALSCHEME_BLFM_TORPOL || defined QUICC_SPATIALSCHEME_WLFL_TORPOL || defined QUICC_SPATIALSCHEME_WLFM_TORPOL
   #include "Diagnostics/SphericalTorPolWrapper.hpp"
   #include "Diagnostics/SphereCflWrapper.hpp"
#elif defined QUICC_SPATIALSCHEME_FFF || defined QUICC_SPATIALSCHEME_TFF || defined QUICC_SPATIALSCHEME_TFT || defined QUICC_SPATIALSCHEME_TTT 
   #include "Diagnostics/StreamVerticalWrapper.hpp"
   #include "Diagnostics/CartesianCflWrapper.hpp"
#endif //defined QUICC_SPATIALSCHEME_TFF_TORPOL 

namespace QuICC {

namespace Diagnostics {

   DiagnosticCoordinator::DiagnosticCoordinator()
      : mcMaxStep(0.1), mcMinStep(1e-10), mFixedStep(-1), mMaxError(-1.0), mCfl(0.0), mStartTime(0.0), mStartTimestep(0.0)
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

      #if defined QUICC_SPATIALSCHEME_TFF_TORPOL
      // Create a cartesian toroidal/poloidal warpper
      } else if(vectors.count(PhysicalNames::VELOCITY))
      {
         SharedCartesianTorPolWrapper spVelocity = SharedCartesianTorPolWrapper(new CartesianTorPolWrapper(vectors.find(PhysicalNames::VELOCITY)->second));

         this->mspCflWrapper = SharedCartesianCflWrapper(new CartesianCflWrapper(spVelocity));

         this->mFixedStep = tstep(1);

      #elif defined QUICC_SPATIALSCHEME_SLFL_TORPOL || defined QUICC_SPATIALSCHEME_SLFM_TORPOL
      // Create a toroidal/poloidal spherical shell wrapper
      } else if(vectors.count(PhysicalNames::VELOCITY))
      {
         SharedSphericalTorPolWrapper spVelocity = SharedSphericalTorPolWrapper(new SphericalTorPolWrapper(vectors.find(PhysicalNames::VELOCITY)->second));

         this->mspCflWrapper = SharedShellCflWrapper(new ShellCflWrapper(spVelocity));

         this->mFixedStep = tstep(1);

      #elif defined QUICC_SPATIALSCHEME_BLFL_TORPOL || defined QUICC_SPATIALSCHEME_BLFM_TORPOL || defined QUICC_SPATIALSCHEME_WLFL_TORPOL || defined QUICC_SPATIALSCHEME_WLFM_TORPOL
      // Create a full sphere wrapper
      } else if(vectors.count(PhysicalNames::VELOCITY))
      {
         SharedSphericalTorPolWrapper spVelocity = SharedSphericalTorPolWrapper(new SphericalTorPolWrapper(vectors.find(PhysicalNames::VELOCITY)->second));

         this->mspCflWrapper = SharedSphereCflWrapper(new SphereCflWrapper(spVelocity));

         this->mFixedStep = tstep(1);

      #elif defined QUICC_SPATIALSCHEME_FFF || defined QUICC_SPATIALSCHEME_TFF || defined QUICC_SPATIALSCHEME_TFT || defined QUICC_SPATIALSCHEME_TTT 
      // Create a stream function and vertical velocity wrapper
      } else if(scalars.count(PhysicalNames::STREAMFUNCTION) && scalars.count(PhysicalNames::VELOCITYZ))
      {
         SharedStreamVerticalWrapper spVelocity = SharedStreamVerticalWrapper(new StreamVerticalWrapper(scalars.find(PhysicalNames::STREAMFUNCTION)->second, scalars.find(PhysicalNames::VELOCITYZ)->second));

         this->mspCflWrapper = SharedCartesianCflWrapper(new CartesianCflWrapper(spVelocity));

         this->mFixedStep = tstep(1);

      #endif //defined QUICC_SPATIALSCHEME_TFF_TORPOL 

      // Required wrapper is not implemented
      } else
      {
         this->mFixedStep = tstep(1);
      }

      if(this->mspCflWrapper)
      {
         this->mspCflWrapper->init(mesh);
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
      #ifdef QUICC_MPI

      if(this->mFixedStep <= 0 && this->mspCflWrapper)
      {
         // Reduce CFL on all CPUs to the global minimum
         MPI_Allreduce(MPI_IN_PLACE, &this->mCfl, 1, Parallel::MpiTypes::type<MHDFloat>(), MPI_MIN, MPI_COMM_WORLD);
      }

      //
      // End of MPI block
      //
      #endif // QUICC_MPI
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
