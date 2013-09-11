/** 
 * @file MpiFramework.cpp
 * @brief Source of the implementation of an MPI framework
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 * @version 0.9.0
 * @date 2013-09-11
 */

// Configuration includes
//
#include "Profiler/ProfilerMacro.h"

// System includes
//
#include <mpi.h>

// External includes
//

// Class include
//
#include "Framework/MpiFramework.hpp"

// Project includes
//
#include "Base/Precision.hpp"

namespace GeoMHDiSCC {

   void MpiFramework::init()
   {
      // Initialise the profiler if needed
      ProfilerMacro_init();

      // Initialise the precision framework
      Precision::init();

      // Initialise MPI system
      MPI_Init(0, 0);

      // Set MPI rank of local CPU
      MPI_Comm_rank(MPI_COMM_WORLD, &MpiFramework::mCpuId);
   }

   void MpiFramework::setup(const int nCpu)
   {
      // Set the number of CPUs
      MpiFramework::mNCpu = nCpu;

      // Get MPI size
      int size;
      MPI_Comm_size(MPI_COMM_WORLD, &size);

      // Check that the framework setup is right
      MpiFramework::checkFramework(size);
   }

   void MpiFramework::synchronize()
   {
      // Create MPI barrier to force synchronisation
      MPI_Barrier(MPI_COMM_WORLD);
   }

   void MpiFramework::finalize()
   {
      // Make sure all finished and are synchronised
      MPI_Barrier(MPI_COMM_WORLD);

      // Finalise MPI system
      MPI_Finalize();
   }

}
