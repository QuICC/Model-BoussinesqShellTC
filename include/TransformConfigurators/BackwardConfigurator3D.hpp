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
#include "TransformConfigurators/TransformTree.hpp"
#include "TransformConfigurators/BackwardConfigurator2D.hpp"
#include "TransformConfigurators/TransformStepsMacro.h"

namespace QuICC {

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
          */
         static void project2D(const TransformTreeEdge& edge, TransformCoordinatorType& coord);

         /**
          * @brief Empty constructor
          */
         BackwardConfigurator3D() {};

         /**
          * @brief Empty destructor
          */
         virtual ~BackwardConfigurator3D() {};

      private: 
   };

#ifdef QUICC_SPATIALDIMENSION_3D
   typedef BackwardConfigurator3D BackwardConfigurator;
#endif //QUICC_SPATIALDIMENSION_3D

}
}

#endif // BACKWARDCONFIGURATOR3D_HPP
