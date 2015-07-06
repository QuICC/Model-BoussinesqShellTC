/** 
 * @file IntegratorTree3DTools.hpp
 * @brief Implementation of tools to work with the integrator transform trees for 3D space
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef INTEGRATORTREE3DTOOLS_HPP
#define INTEGRATORTREE3DTOOLS_HPP

// System includes
//
#include <vector>
#include <map>

// External includes
//

// Project includes
//
#include "TransformConfigurators/IntegratorBranch3D.hpp"
#include "TransformConfigurators/IntegratorTree3D.hpp"

namespace GeoMHDiSCC {

namespace Transform {

   /**
    * @brief Implementation of requirement tools for the integrator transform tree for 3D space
    */
   class IntegratorTree3DTools
   {
      public:
         /**
          * @brief Create vector of projector trees from map of branches
          *
          * @param rTrees     Output integrator trees
          * @param branches   Input tree branches
          */
         static void generateTrees(std::vector<IntegratorTree3D>& rTrees, const std::map<PhysicalNames::Id, std::vector<IntegratorBranch3D> >& branches);
         
      protected:

      private:
         /**
          * @brief Constructor
          */
         IntegratorTree3DTools();

         /**
          * @brief Destructor
          */
         ~IntegratorTree3DTools();
   };

}
}

#endif // INTEGRATORTREE3DTOOLS_HPP
