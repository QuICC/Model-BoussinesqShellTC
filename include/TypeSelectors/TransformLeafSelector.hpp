/** 
 * @file TransformLeafSelector.hpp
 * @brief Typedefs to setup the transform tree leaves
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef TRANSFORMLEAFSELECTOR_HPP
#define TRANSFORMLEAFSELECTOR_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "TypeSelectors/TransformSelector.hpp"
#include "TransformConfigurators/TransformEdge.hpp"

namespace GeoMHDiSCC {

   namespace Transform {

      #if defined GEOMHDISCC_SPATIALSCHEME_TF || defined GEOMHDISCC_SPATIALSCHEME_TT || defined GEOMHDISCC_SPATIALSCHEME_AF || defined GEOMHDISCC_SPATIALSCHEME_CF
         /// Typedefs to simplify definition of integration operators
         typedef TransformSelector<Dimensions::Transform::TRA1D>::Type::IntegratorType IntgSpecType;
         typedef TransformSelector<Dimensions::Transform::TRA2D>::Type::IntegratorType IntgphysType;

         typedef TransformSelector<Dimensions::Transform::TRA1D>::Type::IntegratorType::Id IntgSpecId;
         typedef TransformSelector<Dimensions::Transform::TRA2D>::Type::IntegratorType::Id IntgphysId;

         typedef TransformEdge<IntgPhysId,IntgSpecId,void> IntegratorPhysEdge;
         typedef TransformEdge<IntgSpecId,void,void> IntegratorSpecEdge;

         typedef std::vector<IntegratorSpecEdge>::const_iterator  IntegratorSpecEdge_iterator;
         typedef std::vector<IntegratorPhysEdge>::const_iterator  IntegratorPhysEdge_iterator;

         typedef std::pair<IntegratorSpecEdge_iterator,IntegratorSpecEdge_iterator> IntegratorSpecEdge_range;
         typedef std::pair<IntegratorPhysEdge_iterator,IntegratorPhysEdge_iterator> IntegratorPhysEdge_range;

         /// Typedefs to simplify definition of projection operators
         typedef TransformSelector<Dimensions::Transform::TRA1D>::Type::ProjectorType ProjSpecType;
         typedef TransformSelector<Dimensions::Transform::TRA2D>::Type::ProjectorType ProjPhysType;

         typedef TransformSelector<Dimensions::Transform::TRA1D>::Type::ProjectorType::Id ProjSpecId;
         typedef TransformSelector<Dimensions::Transform::TRA2D>::Type::ProjectorType::Id ProjPhysId;

         typedef TransformEdge<ProjSpecId,ProjPhysId,void> ProjectorSpecEdge;
         typedef TransformEdge<ProjPhysId,void,,void> ProjectorPhysEdge;

         typedef std::vector<ProjectorSpecEdge>::const_iterator  ProjectorSpecEdge_iterator;
         typedef std::vector<ProjectorPhysEdge>::const_iterator  ProjectorPhysEdge_iterator;

         typedef std::pair<ProjectorSpecEdge_iterator,ProjectorSpecEdge_iterator> ProjectorSpecEdge_range;
         typedef std::pair<ProjectorPhysEdge_iterator,ProjectorPhysEdge_iterator> ProjectorPhysEdge_range;
      #else
         /// Typedefs to simplify definition of integration operators
         typedef TransformSelector<Dimensions::Transform::TRA1D>::Type::IntegratorType IntgSpecType;
         typedef TransformSelector<Dimensions::Transform::TRA2D>::Type::IntegratorType IntgPartType;
         typedef TransformSelector<Dimensions::Transform::TRA3D>::Type::IntegratorType IntgPhysType;

         typedef TransformSelector<Dimensions::Transform::TRA1D>::Type::IntegratorType::Id IntgSpecId;
         typedef TransformSelector<Dimensions::Transform::TRA2D>::Type::IntegratorType::Id IntgPartId;
         typedef TransformSelector<Dimensions::Transform::TRA3D>::Type::IntegratorType::Id IntgPhysId;

         typedef TransformEdge<IntgPhysId,IntgPartId,IntgSpecId> IntegratorPhysEdge;
         typedef TransformEdge<IntgPartId,IntgSpecId,void> IntegratorPartEdge;
         typedef TransformEdge<IntgSpecId,void,void> IntegratorSpecEdge;

         typedef std::vector<IntegratorSpecEdge>::const_iterator  IntegratorSpecEdge_iterator;
         typedef std::vector<IntegratorPartEdge>::const_iterator  IntegratorPartEdge_iterator;
         typedef std::vector<IntegratorPhysEdge>::const_iterator  IntegratorPhysEdge_iterator;

         typedef std::pair<IntegratorSpecEdge_iterator,IntegratorSpecEdge_iterator> IntegratorSpecEdge_range;
         typedef std::pair<IntegratorPartEdge_iterator,IntegratorPartEdge_iterator> IntegratorPartEdge_range;
         typedef std::pair<IntegratorPhysEdge_iterator,IntegratorPhysEdge_iterator> IntegratorPhysEdge_range;

         /// Typedefs to simplify definition of projection operators
         typedef TransformSelector<Dimensions::Transform::TRA1D>::Type::ProjectorType ProjSpecType;
         typedef TransformSelector<Dimensions::Transform::TRA2D>::Type::ProjectorType ProjPartType;
         typedef TransformSelector<Dimensions::Transform::TRA3D>::Type::ProjectorType ProjPhysType;

         typedef TransformSelector<Dimensions::Transform::TRA1D>::Type::ProjectorType::Id ProjSpecId;
         typedef TransformSelector<Dimensions::Transform::TRA2D>::Type::ProjectorType::Id ProjPartId;
         typedef TransformSelector<Dimensions::Transform::TRA3D>::Type::ProjectorType::Id ProjPhysId;

         typedef TransformEdge<ProjSpecId,ProjPartId,ProjPhysId> ProjectorSpecEdge;
         typedef TransformEdge<ProjPartId,ProjPhysId,void> ProjectorPartEdge;
         typedef TransformEdge<ProjPhysId,void,void> ProjectorPhysEdge;

         typedef std::vector<ProjectorSpecEdge>::const_iterator  ProjectorSpecEdge_iterator;
         typedef std::vector<ProjectorPartEdge>::const_iterator  ProjectorPartEdge_iterator;
         typedef std::vector<ProjectorPhysEdge>::const_iterator  ProjectorPhysEdge_iterator;

         typedef std::pair<ProjectorSpecEdge_iterator,ProjectorSpecEdge_iterator> ProjectorSpecEdge_range;
         typedef std::pair<ProjectorPartEdge_iterator,ProjectorPartEdge_iterator> ProjectorPartEdge_range;
         typedef std::pair<ProjectorPhysEdge_iterator,ProjectorPhysEdge_iterator> ProjectorPhysEdge_range;
      #endif //defined GEOMHDISCC_SPATIALSCHEME_TF || defined GEOMHDISCC_SPATIALSCHEME_TF || defined GEOMHDISCC_SPATIALSCHEME_AF || defined GEOMHDISCC_SPATIALSCHEME_CF
   }
}

#endif // TRANSFORMLEAFSELECTOR_HPP
