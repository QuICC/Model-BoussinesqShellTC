/** \file MpiFramework.hpp
 *  \brief Implementation of the MPI framework
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
   };

}

#endif // MPIFRAMEWORK_HPP
