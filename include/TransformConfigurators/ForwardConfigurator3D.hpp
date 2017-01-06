/** 
 * @file ForwardConfigurator3D.hpp
 * @brief This class defines the base operations for a forward transform in 3D space
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef FORWARDCONFIGURATOR3D_HPP
#define FORWARDCONFIGURATOR3D_HPP

// Configuration includes
// 
#include "Profiler/ProfilerMacro.h"

// System includes
//

// External includes
//

// Project includes
//
#include "Equations/IScalarEquation.hpp"
#include "Equations/IVectorEquation.hpp"
#include "TypeSelectors/VariableSelector.hpp"
#include "TypeSelectors/TransformCommSelector.hpp"
#include "TransformConfigurators/TransformTree.hpp"
#include "TransformConfigurators/ForwardConfigurator2D.hpp"
#include "TransformConfigurators/TransformStepsMacro.h"

namespace GeoMHDiSCC {

namespace Transform {

   /**
    * @brief This class defines the base operations for a forward transform in 3D space
    */
   class ForwardConfigurator3D: public ForwardConfigurator2D
   {
      public:

      protected:
         /**
          * @brief Compute the integration transform of the second dimension
          *
          * @param coord   Transform coordinator
          */
         static void integrate2D(const TransformTreeEdge& edge, TransformCoordinatorType& coord);

         /**
          * @brief Empty constructor
          */
         ForwardConfigurator3D() {};

         /**
          * @brief Empty destructor
          */
         virtual ~ForwardConfigurator3D() {};

      private: 
   };

#ifdef QUICC_SPATIALDIMENSION_3D
   typedef ForwardConfigurator3D ForwardConfigurator;
#endif //QUICC_SPATIALDIMENSION_3D

}
}

#endif // FORWARDCONFIGURATOR3D_HPP
