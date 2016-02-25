/** 
 * @file IntegratorTree2DTools.hpp
 * @brief Implementation of tools to work with the integrator transform trees for 2D space
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef INTEGRATORTREE2DTOOLS_HPP
#define INTEGRATORTREE2DTOOLS_HPP

// System includes
//
#include <vector>
#include <map>

// External includes
//

// Project includes
//
#include "TransformConfigurators/IntegratorBranch2D.hpp"
#include "TransformConfigurators/IntegratorTree2D.hpp"

namespace GeoMHDiSCC {

namespace Transform {

   /**
    * @brief Implementation of requirement tools for the integrator transform tree for 2D space
    */
   class IntegratorTree2DTools
   {
      public:
         /**
          * @brief Create vector of projector trees from map of branches
          *
          * @param rTrees     Output integrator trees
          * @param branches   Input tree branches
          */
         static void generateTrees(std::vector<IntegratorTree2D>& rTrees, const std::map<PhysicalNames::Id, std::vector<IntegratorBranch2D> >& branches);
         
      protected:

      private:
         /**
          * @brief Constructor
          */
         IntegratorTree2DTools();

         /**
          * @brief Destructor
          */
         ~IntegratorTree2DTools();
   };

}
}

#endif // INTEGRATORTREE2DTOOLS_HPP
