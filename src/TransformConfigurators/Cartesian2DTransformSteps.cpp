/** 
 * @file Cartesian2DTransformSteps.cpp
 * @brief Source of the implementation of the 2D cartesian transform steps
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "TransformConfigurators/Cartesian2DTransformSteps.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Transform {

namespace TransformSteps {

   std::vector<IntegratorBranch2D>  forwardScalar(const std::vector<std::pair<FieldComponents::Spectral::Id,int> >& components, const bool isNL)
   {
      assert(components.size() == 1);
      std::vector<IntegratorBranch2D> transform;

      FieldComponents::Spectral::Id scalId = components.at(0).first;

      if(isNL)
      {
         transform.push_back(IntegratorBranch2D(FieldComponents::Physical::SCALAR, IntgPhysType::INTG, IntgSpecType::INTG, scalId, FieldType::SCALAR));
      } else
      {
         transform.push_back(IntegratorBranch2D(FieldComponents::Physical::SCALAR, IntgPhysType::INTG, IntgSpecType::INTG, scalId, FieldType::SCALAR));
      }

      return transform;
   }

   std::vector<IntegratorBranch2D>  forwardVector(const std::vector<std::pair<FieldComponents::Spectral::Id,int> >& components, const bool isNL)
   {
      assert(components.size() == 2);
      std::vector<IntegratorBranch2D> transform;

      if(isNL)
      {
         transform.push_back(IntegratorBranch2D(FieldComponents::Physical::ONE, IntgPhysType::INTG, IntgSpecType::INTG, FieldComponents::Spectral::ONE, FieldType::VECTOR));

         transform.push_back(IntegratorBranch2D(FieldComponents::Physical::TWO, IntgPhysType::INTG, IntgSpecType::INTG, FieldComponents::Spectral::TWO, FieldType::VECTOR));
      } else
      {
         transform.push_back(IntegratorBranch2D(FieldComponents::Physical::ONE, IntgPhysType::INTG, IntgSpecType::INTG, FieldComponents::Spectral::ONE, FieldType::VECTOR));

         transform.push_back(IntegratorBranch2D(FieldComponents::Physical::TWO, IntgPhysType::INTG, IntgSpecType::INTG, FieldComponents::Spectral::TWO, FieldType::VECTOR));
      }

      return transform;
   }

   std::vector<ProjectorBranch2D>  backwardScalar(const std::map<FieldComponents::Physical::Id,bool>& req)
   {
      std::vector<ProjectorBranch2D> transform;

      transform.push_back(ProjectorBranch2D(FieldComponents::Spectral::SCALAR, ProjSpecType::PROJ, ProjPhysType::PROJ, FieldComponents::Physical::SCALAR, FieldType::SCALAR));

      return transform;
   }

   std::vector<ProjectorBranch2D>  backwardGradient(const std::map<FieldComponents::Physical::Id,bool>& req)
   {
      std::vector<ProjectorBranch2D> transform;

      if(req.find(FieldComponents::Physical::ONE)->second)
      {
         transform.push_back(ProjectorBranch2D(FieldComponents::Spectral::SCALAR, ProjSpecType::DIFF, ProjPhysType::PROJ, FieldComponents::Physical::ONE, FieldType::GRADIENT));
      }

      if(req.find(FieldComponents::Physical::TWO)->second)
      {
         transform.push_back(ProjectorBranch2D(FieldComponents::Spectral::SCALAR, ProjSpecType::PROJ, ProjPhysType::DIFF, FieldComponents::Physical::TWO, FieldType::GRADIENT));
      }

      return transform;
   }

   std::vector<ProjectorBranch2D>  backwardGradient2(FieldComponents::Spectral::Id id, const std::map<std::pair<FieldComponents::Physical::Id,FieldComponents::Physical::Id>,bool>& req)
   {
      std::vector<ProjectorBranch2D> transform;
      std::pair<FieldComponents::Physical::Id,FieldComponents::Physical::Id>  pairId;

      pairId = std::make_pair(FieldComponents::Physical::ONE,FieldComponents::Physical::ONE);
      if(req.find(pairId)->second)
      {
         transform.push_back(ProjectorBranch2D(FieldComponents::Spectral::SCALAR, ProjSpecType::DIFF2, ProjPhysType::PROJ, pairId, FieldType::GRADIENT2));
      }

      pairId = std::make_pair(FieldComponents::Physical::ONE,FieldComponents::Physical::TWO);
      if(req.find(pairId)->second)
      {
         transform.push_back(ProjectorBranch2D(FieldComponents::Spectral::SCALAR, ProjSpecType::DIFF, ProjPhysType::DIFF, pairId, FieldType::GRADIENT2));
      }

      pairId = std::make_pair(FieldComponents::Physical::TWO,FieldComponents::Physical::TWO);
      if(req.find(pairId)->second)
      {
         transform.push_back(ProjectorBranch2D(FieldComponents::Spectral::SCALAR, ProjSpecType::PROJ, ProjPhysType::DIFF2, pairId, FieldType::GRADIENT2));
      }

      return transform;
   }

   std::vector<ProjectorBranch2D>  backwardVector(const std::map<FieldComponents::Physical::Id,bool>& req)
   {
      std::vector<ProjectorBranch2D> transform;

      if(req.find(FieldComponents::Physical::ONE)->second)
      {
         transform.push_back(ProjectorBranch2D(FieldComponents::Spectral::ONE, ProjSpecType::PROJ, ProjPhysType::PROJ, FieldComponents::Physical::ONE, FieldType::VECTOR));
      }

      if(req.find(FieldComponents::Physical::TWO)->second)
      {
         transform.push_back(ProjectorBranch2D(FieldComponents::Spectral::TWO, ProjSpecType::PROJ, ProjPhysType::PROJ, FieldComponents::Physical::TWO, FieldType::VECTOR));
      }

      return transform;
   }

   std::vector<ProjectorBranch2D>  backwardVGradient(FieldComponents::Spectral::Id id, const std::map<FieldComponents::Physical::Id,bool>& req)
   {
      std::vector<ProjectorBranch2D> transform;

      if(req.find(FieldComponents::Physical::ONE)->second)
      {
         transform.push_back(ProjectorBranch2D(id, ProjSpecType::DIFF, ProjPhysType::PROJ, FieldComponents::Physical::ONE, FieldType::GRADIENT));
      }

      if(req.find(FieldComponents::Physical::TWO)->second)
      {
         transform.push_back(ProjectorBranch2D(id, ProjSpecType::PROJ, ProjPhysType::DIFF, FieldComponents::Physical::TWO, FieldType::GRADIENT));
      }

      return transform;
   }

   std::vector<ProjectorBranch2D>  backwardCurl(const std::map<FieldComponents::Physical::Id,bool>& req)
   {
      throw Exception("Curl cannot be computed in 2D!");

      std::vector<ProjectorBranch2D> transform;

      return transform;
   }

   std::vector<ProjectorBranch2D>  backwardDivergence()
   {
      std::vector<ProjectorBranch2D> transform;

      transform.push_back(ProjectorBranch2D(FieldComponents::Spectral::ONE, ProjSpecType::DIFF, ProjPhysType::PROJ, FieldComponents::Physical::SCALAR, FieldType::DIVERGENCE));

      transform.push_back(ProjectorBranch2D(FieldComponents::Spectral::TWO, ProjSpecType::PROJ, ProjPhysType::DIFF, FieldComponents::Physical::SCALAR, FieldType::DIVERGENCE));

      return transform;
   }

}
}
}
