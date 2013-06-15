/** \file StorageProfiler.hpp
 *  \brief Implementation of a storage profiler
 *
 *  \mhdTodo Cleanup the StorageProfiler class similarly to the Profiler class
 */

#ifndef STORAGEPROFILER_HPP
#define STORAGEPROFILER_HPP

// Configuration includes
//

// System includes
//
#include <map>

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"

namespace GeoMHDiSCC {

namespace Debug {

   /**
    * \brief Implementation of a storage profiler
    */
   class StorageProfiler
   {
      public:
         /**
          * @brief List of storage locations
          */
         enum StoragePoint {VARIABLES, TRANSFORMS, TEMPORARIES, DIAGNOSTICEQUATION, PROGNOSTICEQUATION, MPI, IO, SCALAR1, SCALAR2, SCALAR3, VECTOR1, VECTOR2, VECTOR3, TRANSFORM1D, TRANSFORM2D, TRANSFORM3D, TRANSFORM1DTEMP, TRANSFORM2DTEMP, TRANSFORM3DTEMP, MPIBUFFERS, MPITYPES, MPICOMM, NMAXSTORAGEPOINT};

         /**
          * @brief Update required storage
          *
          * @param point   Storage location
          * @param memory  Required storage
          */
         static void update(StoragePoint point, MHDFloat memory);

         /**
          * @brief Analyze the measured storage requirements
          *
          * @param min  Array of MPI minimal values
          * @param max  Array of MPI maximal values
          */
         static void analyze(Array& min, Array& max);

         /**
          * @brief Create enum ID to name string map
          *
          * @param map  Storage for the map to build
          */
         static void createNameMap(std::map<StoragePoint, std::string>& map);

         /**
          * @brief Storage requirements
          */
         static std::map<StoragePoint, MHDFloat> requirements;
         
      protected:

      private:
         /**
          * @brief Constructor
          */
         StorageProfiler();

         /**
          * @brief Destructor
          */
         ~StorageProfiler();
   };

}
}

#endif // STORAGEPROFILER_HPP
