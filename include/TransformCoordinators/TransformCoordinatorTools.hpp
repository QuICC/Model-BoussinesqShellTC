/** 
 * @file TransformCoordinatorTools.hpp
 * @brief Implementation of transform coordinator tools
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef TRANSFORMCOORDINATORTOOLS_HPP
#define TRANSFORMCOORDINATORTOOLS_HPP

// System includes
//
#include <map>

// External includes
//

// Project includes
//
#include "TypeSelectors/TransformSelector.hpp"
#include "TransformGroupers/IForwardGrouper.hpp"
#include "TransformGroupers/IBackwardGrouper.hpp"
#include "TransformConfigurators/ProjectorTree.hpp"

namespace GeoMHDiSCC {

namespace Transform {

   /**
    * @brief Implementation of transform coordinator tools
    */
   class TransformCoordinatorTools
   {
      public:
         /**
          * @brief Initialise the transform coordinator
          *
          * @param rCoord           Transform coordinator
          * @param spFwdGrouper     Forward transform communication grouper
          * @param spBwdGrouper     Backward transform communication grouper
          * @param integratorTree   Transform integrator tree
          * @param projectorTree    Transform projector tree
          * @param spRes            Shared resolution
          * @param runOptions       Available run options map
          */
         static void init(TransformCoordinatorType& rCoord, SharedIForwardGrouper spFwdGrouper, SharedIBackwardGrouper spBwdGrouper, const std::vector<Transform::IntegratorTree>& integratorTree, const std::vector<Transform::ProjectorTree>& projectorTree, SharedResolution spRes, const std::map<NonDimensional::Id,MHDFloat>& runOptions);
         
      protected:

      private:
         /**
          * @brief Constructor
          */
         TransformCoordinatorTools();

         /**
          * @brief Destructor
          */
         ~TransformCoordinatorTools();
   };

}
}

#endif // TRANSFORMCOORDINATORTOOLS_HPP
