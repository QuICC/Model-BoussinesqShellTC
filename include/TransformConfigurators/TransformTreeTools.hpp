/** 
 * @file TransformTreeTools.hpp
 * @brief Implementation of tools to work with the transform trees
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef TRANSFORMTREETOOLS_HPP
#define TRANSFORMTREETOOLS_HPP

// System includes
//
#include <vector>
#include <map>
#include <string>

// External includes
//

// Project includes
//
#include "Enums/TransformDirection.hpp"
#include "TransformConfigurators/TransformPath.hpp"
#include "TransformConfigurators/TransformTree.hpp"

namespace GeoMHDiSCC {

namespace Transform {

   /**
    * @brief Implementation of tools to work with the projector transform trees for 3D space
    */
   class TransformTreeTools
   {
      public:
         /**
          * @brief Create vector of transform trees from vector of transform paths
          *
          * @param rTrees     Output transform trees
          * @param branches   Input transfrom paths
          * @param dir        Direction of transform (only needed for graph output)
          */
         static void generateTrees(std::vector<TransformTree>& rTrees, const std::map<PhysicalNames::Id, std::vector<TransformPath> >& branches, const TransformDirection::Id dir, const std::string& prepend = "");
         
      protected:

      private:
         typedef std::pair<FieldType::Id,std::vector<int> >  OutKey;
         typedef std::set<OutKey> SetOutIds;
         typedef std::tr1::tuple<int,FieldType::Id,std::vector<int> > OptKey;
         typedef std::map<OptKey,std::tr1::tuple<int,int,bool> > MapOptIds;

         /**
          * @brief Recursive grows of the tree
          */
         static void growTree(std::map<PhysicalNames::Id, std::vector<TransformPath> >::const_iterator nameIt, std::set<int>::const_iterator compIt, const int dim, std::vector<int>& path, TransformTreeEdge& rPrev);

         /**
          * @brief Build the trees recursively
          */
         static void buildTrees(std::vector<TransformTree>& rTrees, const std::map<PhysicalNames::Id, std::vector<TransformPath> >& branches);

         /**
          * @brief Set edge input information
          */
         static void setInputInfoEdge(TransformTreeEdge& edge);

         /**
          * @brief Set edge output order
          */
         static void setOutputOrder(TransformTreeEdge& edge, SetOutIds& outIds);

         /**
          * @brief Optimize output count
          */
         static void optimizeOutputCount(TransformTreeEdge& edge, MapOptIds& optIds, SetOutIds& outIds, int& counter);

         /**
          * @brief Optimize output by pruning extra outputs
          */
         static void pruneOutput(TransformTreeEdge& edge, MapOptIds& optIds);

         /**
          * @brief Finalize the trees
          */
         static void finalizeTrees(std::vector<TransformTree>& rTrees);

         /**
          * @brief Finalize the trees
          */
         static void optimizeTrees(std::vector<TransformTree>& rTrees);

         /**
          * @brief Constructor
          */
         TransformTreeTools();

         /**
          * @brief Destructor
          */
         ~TransformTreeTools();
   };

}
}

#endif // TRANSFORMTREETOOLS_HPP
