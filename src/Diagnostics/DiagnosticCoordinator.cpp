/**
 * @file DiagnosticCoordinator.cpp
 * @brief Source of the diagnostic coordinator
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Configuration includes
//
#include "Framework/FrameworkMacro.h"

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
      : MAXSTEP_LOCATION(-101), MINSTEP_LOCATION(-102), FIXEDSTEP_LOCATION(-100), mcMaxStep(0.1), mcMinStep(1e-10), mFixedStep(-1), mMaxError(-1.0), mCfl(2,1), mStartTime(0.0), mStartTimestep(0.0)
   {
      this->mCfl.setZero();
   }

   DiagnosticCoordinator::~DiagnosticCoordinator()
   {
   }

   void DiagnosticCoordinator::init(const std::vector<Array>& mesh, const std::map<PhysicalNames::Id, Datatypes::SharedScalarVariableType>&  scalars, const std::map<PhysicalNames::Id, Datatypes::SharedVectorVariableType>&  vectors, const Array& tstep, const std::map<NonDimensional::Id,MHDFloat>& params)
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
      } else if(vectors.count(PhysicalNames::VELOCITY) && vectors.count(PhysicalNames::MAGNETIC))
      {
         SharedSphericalTorPolWrapper spVelocity = SharedSphericalTorPolWrapper(new SphericalTorPolWrapper(vectors.find(PhysicalNames::VELOCITY)->second));
         SharedSphericalTorPolWrapper spMagnetic = SharedSphericalTorPolWrapper(new SphericalTorPolWrapper(vectors.find(PhysicalNames::MAGNETIC)->second));

         this->mspCflWrapper = SharedShellCflWrapper(new ShellCflWrapper(spVelocity, spMagnetic, params));

         this->mFixedStep = tstep(1);

      // Create a toroidal/poloidal spherical shell wrapper
      } else if(vectors.count(PhysicalNames::VELOCITY))
      {
         SharedSphericalTorPolWrapper spVelocity = SharedSphericalTorPolWrapper(new SphericalTorPolWrapper(vectors.find(PhysicalNames::VELOCITY)->second));

         this->mspCflWrapper = SharedShellCflWrapper(new ShellCflWrapper(spVelocity, params));

         this->mFixedStep = tstep(1);

      #elif defined QUICC_SPATIALSCHEME_BLFL_TORPOL || defined QUICC_SPATIALSCHEME_BLFM_TORPOL || defined QUICC_SPATIALSCHEME_WLFL_TORPOL || defined QUICC_SPATIALSCHEME_WLFM_TORPOL
      // Create a full sphere wrapper
      } else if(vectors.count(PhysicalNames::VELOCITY) && vectors.count(PhysicalNames::MAGNETIC))
      {
         SharedSphericalTorPolWrapper spVelocity = SharedSphericalTorPolWrapper(new SphericalTorPolWrapper(vectors.find(PhysicalNames::VELOCITY)->second));
         SharedSphericalTorPolWrapper spMagnetic = SharedSphericalTorPolWrapper(new SphericalTorPolWrapper(vectors.find(PhysicalNames::MAGNETIC)->second));

         this->mspCflWrapper = SharedSphereCflWrapper(new SphereCflWrapper(spVelocity, spMagnetic, params));

         this->mFixedStep = tstep(1);
      // Create a full sphere wrapper
      } else if(vectors.count(PhysicalNames::VELOCITY))
      {
         SharedSphericalTorPolWrapper spVelocity = SharedSphericalTorPolWrapper(new SphericalTorPolWrapper(vectors.find(PhysicalNames::VELOCITY)->second));

         this->mspCflWrapper = SharedSphereCflWrapper(new SphereCflWrapper(spVelocity, params));

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
         this->mCfl(0,0) = this->mcMinStep;
         this->mCfl(1,0) = MINSTEP_LOCATION;

      } else if(this->mFixedStep > 0)
      {
         this->mCfl(0,0) = this->mFixedStep;
         this->mCfl(1,0) = FIXEDSTEP_LOCATION;

      // Compute initial CFL condition
      } else if(this->mspCflWrapper)
      {
         // Compute CFL for initial state
         this->mCfl = this->mspCflWrapper->initialCfl();

         if(this->mcMinStep < this->mCfl(0,0))
         {
            this->mCfl(0,0) = this->mcMinStep;
            this->mCfl(1,0) = MINSTEP_LOCATION;
         }
      }
   }

   void DiagnosticCoordinator::updateCfl()
   {
      // Used fixed timestep
      if(this->mFixedStep > 0)
      {
         this->mCfl(0,0) = this->mFixedStep;
         this->mCfl(1,0) = FIXEDSTEP_LOCATION;

      // Compute CFL condition
      } else if(this->mspCflWrapper)
      {
         // Safety assert
         assert(this->mspCflWrapper);

         this->mCfl = this->mspCflWrapper->cfl();

         // Check for maximum timestep
         if(this->mcMaxStep < this->mCfl(0,0))
         {
            this->mCfl(0,0) = this->mcMaxStep;
            this->mCfl(1,0) = MAXSTEP_LOCATION;
         }

         if(-this->mFixedStep < this->mCfl(0,0))
         {
            this->mCfl(0,0) = -this->mFixedStep;
            this->mCfl(1,0) = FIXEDSTEP_LOCATION;
         }
      }
   }

   void DiagnosticCoordinator::synchronize()
   {
      // Start of MPI block
      #ifdef QUICC_MPI

      if(this->mFixedStep <= 0 && this->mspCflWrapper)
      {
         // Create MPI operation
         MPI_Op op;
         MPI_Op_create(mpi_cfl_min, true, &op);

         // Create MPI datatype
         MPI_Datatype ctype;
         MPI_Type_contiguous(2, Parallel::MpiTypes::type<MHDFloat>(), &ctype);
         MPI_Type_commit(&ctype);

         // Reduce CFL on all CPUs to the global minimum
         MPI_Allreduce(MPI_IN_PLACE, this->mCfl.data(), this->mCfl.cols(), ctype, op, MPI_COMM_WORLD);
      }

      // End of MPI block
      #endif // QUICC_MPI
   }

   MHDFloat DiagnosticCoordinator::maxError() const
   {
      return this->mMaxError;
   }

   const Matrix& DiagnosticCoordinator::cfl() const
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

   #ifdef QUICC_MPI
   void mpi_cfl_min(void* a, void* b, int* len, MPI_Datatype* type)
   {
      MHDFloat* in = static_cast<MHDFloat*>(a);
      MHDFloat* inout = static_cast<MHDFloat*>(b);

      int cols = *len;
      int rows;
      MPI_Type_size(*type, &rows);
      rows /= sizeof(MHDFloat);

      for (int i = 0; i < cols; i++)
      {
         if(in[i*rows] < inout[i*rows])
         {
            inout[i*rows] = in[i*rows];
            for(int j = 1; j < rows; j++)
            {
               inout[i*rows+j] = in[i*rows+j];
            }
         }
      }
   }
   #endif //QUICC_MPI

}
}
