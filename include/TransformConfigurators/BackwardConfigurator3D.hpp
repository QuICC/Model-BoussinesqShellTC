/** 
 * @file BackwardConfigurator3D.hpp
 * @brief This class defines the base operations for a backward transform in 3D
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef BACKWARDCONFIGURATOR3D_HPP
#define BACKWARDCONFIGURATOR3D_HPP

// Configuration includes
// 
#include "Profiler/ProfilerMacro.h"

// System includes
//

// External includes
//

// Project includes
//
#include "TypeSelectors/VariableSelector.hpp"
#include "TypeSelectors/TransformCommSelector.hpp"
#include "TypeSelectors/TransformTreeSelector.hpp"
#include "TransformConfigurators/BackwardConfigurator2D.hpp"
#include "TransformConfigurators/TransformStepsMacro.h"

namespace GeoMHDiSCC {

namespace Transform {

   /**
    * @brief This class defines the base operations for a backward transform in 3D
    */
   class BackwardConfigurator3D: public BackwardConfigurator2D
   {
      public:

      protected:
         /**
          * @brief Compute the projection transform of the second dimension
          *
          * @param edge    Transform tree edge
          * @param coord   Transform coordinator
          * @param recover Recover input data?
          * @param hold    Hold input data?
          */
         static void project2D(const ProjectorPartEdge& edge, TransformCoordinatorType& coord, const bool recover, const bool hold);

         /**
          * @brief Empty constructor
          */
         BackwardConfigurator3D();

         /**
          * @brief Empty destructor
          */
         virtual ~BackwardConfigurator3D();

      private: 
   };

}
}

#endif // BACKWARDCONFIGURATOR3D_HPP
