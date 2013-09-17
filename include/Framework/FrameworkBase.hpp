/**
 * @file FrameworkBase.hpp
 * @brief Base of the implementation of the framework 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef FRAMEWORKBASE_HPP
#define FRAMEWORKBASE_HPP

// System includes
//
#include <cassert>

// External includes
//

// Project includes
//

namespace GeoMHDiSCC {

   /**
    * @brief This class give basic information required for setting up the framework
    */
   class FrameworkBase
   {
      public:
         /**
          * @brief ID of node allowed to do serial IO
          */
         static const int IO_RANK;

         /**
          * @brief Check if local core is allowed to do IO
          */
         static bool allowsIO();

         /**
          * @brief Get the number of CPUs
          */
         static int nCpu();

         /**
          * @brief Get the local CPU ID
          */
         static int id();
         
      protected:

         /**
          * @brief Number of computation cores
          */
         static int mNCpu;

         /**
          * @brief ID of the current node (ie MPI rank)
          */
         static int mCpuId;

         /**
          * @brief Check setup of framework
          *
          * @param nCpu Number of CPUs that have been requested
          */
         static void checkFramework(const int nCpu);

         /**
          * @brief Constructor
          */
         FrameworkBase();

         /**
          * @brief Destructor
          */
         virtual ~FrameworkBase();

      private:
   };

   inline int FrameworkBase::nCpu()
   {
      // Safety assert to avoid uninitialised use
      assert(FrameworkBase::mNCpu > 0);

      return FrameworkBase::mNCpu;
   }

   inline int FrameworkBase::id()
   {
      // Safety assert to avoid uninitialised use
      assert(FrameworkBase::mCpuId >= 0);

      return FrameworkBase::mCpuId;
   }

}

#endif // FRAMEWORKBASE_HPP
