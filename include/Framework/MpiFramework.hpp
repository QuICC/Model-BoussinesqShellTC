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
#include <vector>
#include <set>
#include <map>

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
          * @brief Synchronise a sub communicator
          */
         static void syncSubComm(const SubCommId id, const int idx);

         /**
          * @brief Get a sub communicator
          */
         static MPI_Comm getSubComm(const SubCommId id, const int idx);

         /**
          * @brief Get a sub group
          */
         static MPI_Group getSubGroup(const SubCommId id, const int idx);

         /**
          * @brief Init the MPI sub group and sub communicator vectors
          */
         static void initSubComm(const SubCommId id, const int size);

         /**
          * @brief Set the MPI sub group and sub communicator
          */
         static void setSubComm(const SubCommId id, const int idx, const std::set<int>& ranks);
         
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
          * @brief Storage for the MPI sub groups
          */
         static std::map<SubCommId, std::vector<MPI_Group> > mSubGroup;

         /**
          * @brief Storage for the spectral MPI communicator
          */
         static std::map<SubCommId, std::vector<MPI_Comm> > mSubComm;
   };

}

#endif // MPIFRAMEWORK_HPP
