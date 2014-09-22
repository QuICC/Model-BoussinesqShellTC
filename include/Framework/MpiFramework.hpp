/**
 * @file MpiFramework.hpp
 * @brief Implementation of the MPI framework 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef MPIFRAMEWORK_HPP
#define MPIFRAMEWORK_HPP

// System includes
//
#include <mpi.h>

// External includes
//

// Project includes
//
#include "Framework/FrameworkBase.hpp"

namespace GeoMHDiSCC {

   /**
    * @brief This class defines the framework for an MPI Parallel code
    */
   class MpiFramework: public FrameworkBase
   {
      public:
         /**
          * @brief Initialise the MPI system
          */
         static void init();

         /**
          * @brief Setup the MPI system
          */
         static void setup(const int nCpu);

         /**
          * @brief Synchronise
          */
         static void synchronize();

         /**
          * @brief Finalise the MPI system
          */
         static void finalize();

         /**
          * @brief Spectral CPUs MPI group
          */
         static MPI_Group spectralGroup();

         /**
          * @brief Spectral CPUs MPI communicator
          */
         static MPI_Comm spectralComm();

         /**
          * @brief Set the spectral MPI group and communicator
          */
         static void setSpectralComm();
         
      protected:

      private:
         /**
          * @brief Constructor
          */
         MpiFramework();

         /**
          * @brief Destructor
          */
         ~MpiFramework();

         /**
          * @brief Storage for the spectral MPI group
          */
         static MPI_Group mSpecGroup;  

         /**
          * @brief Storage for the spectral MPI communicator
          */
         static MPI_Comm mSpecComm;  
   };

}

#endif // MPIFRAMEWORK_HPP
