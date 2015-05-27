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
#include "Base/Typedefs.hpp"

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
          * @brief Add CPU group IDs
          */
         static void addTransformComm(const ArrayI& ids);

         /**
          * @brief Synchronize transform comm
          */
         static void syncTransform(const int traId);

         /**
          * @brief Get transform rank
          */
         static int transformId(const int traId);

         /**
          * @brief Get transform CPU group IDs
          */
         static const ArrayI& transformCpus(const int traId);

         /**
          * @brief Get transform MPI group
          */
         static MPI_Group transformGroup(const int traId);

         /**
          * @brief Get transform MPI comm
          */
         static MPI_Comm transformComm(const int traId);

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
          * @brief Local CPU rank in transform group
          */
         static std::vector<int> mTransformIds;

         /**
          * @brief IDs of the CPUs in transform communication groups
          */
         static std::vector<ArrayI> mTransformCpus;

         /**
          * @brief MPI groups of the CPUs in transform communication groups
          */
         static std::vector<MPI_Group> mTransformGroups;

         /**
          * @brief MPI communicators of the CPUs in transform communication groups
          */
         static std::vector<MPI_Comm> mTransformComms;

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
