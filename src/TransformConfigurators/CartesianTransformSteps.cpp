/** 
 * @file CartesianTransformSteps.cpp
 * @brief Source of the implementation of the cartesian transform steps
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
#include "TransformConfigurators/CartesianTransformSteps.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Transform {

namespace TransformSteps {

   std::vector<IntegratorBranch3D>  forwardScalar(const std::vector<std::pair<FieldComponents::Spectral::Id,int> >& components, const bool isNL)
   {
      assert(components.size() == 1);
      std::vector<IntegratorBranch3D> transform;

      FieldComponents::Spectral::Id scalId = components.at(0).first;

      if(isNL)
      {
         transform.push_back(IntegratorBranch3D(FieldComponents::Physical::SCALAR, IntgPhysType::INTG, IntgPartType::INTG, IntgSpecType::INTG, scalId, FieldType::SCALAR));
      } else
      {
         transform.push_back(IntegratorBranch3D(FieldComponents::Physical::SCALAR, IntgPhysType::INTG, IntgPartType::INTG, IntgSpecType::INTG, scalId, FieldType::SCALAR));
      }

      return transform;
   }

   std::vector<IntegratorBranch3D>  forwardVector(const std::vector<std::pair<FieldComponents::Spectral::Id,int> >& components, const bool isNL)
   {
      assert(components.size() == 3);
      std::vector<IntegratorBranch3D> transform;

      if(isNL)
      {
         transform.push_back(IntegratorBranch3D(FieldComponents::Physical::ONE, IntgPhysType::INTG, IntgPartType::INTG, IntgSpecType::INTG, FieldComponents::Spectral::ONE, FieldType::VECTOR));

         transform.push_back(IntegratorBranch3D(FieldComponents::Physical::TWO, IntgPhysType::INTG, IntgPartType::INTG, IntgSpecType::INTG, FieldComponents::Spectral::TWO, FieldType::VECTOR));

         transform.push_back(IntegratorBranch3D(FieldComponents::Physical::THREE, IntgPhysType::INTG, IntgPartType::INTG, IntgSpecType::INTG, FieldComponents::Spectral::THREE, FieldType::VECTOR));
      } else
      {
         transform.push_back(IntegratorBranch3D(FieldComponents::Physical::ONE, IntgPhysType::INTG, IntgPartType::INTG, IntgSpecType::INTG, FieldComponents::Spectral::ONE, FieldType::VECTOR));

         transform.push_back(IntegratorBranch3D(FieldComponents::Physical::TWO, IntgPhysType::INTG, IntgPartType::INTG, IntgSpecType::INTG, FieldComponents::Spectral::TWO, FieldType::VECTOR));

         transform.push_back(IntegratorBranch3D(FieldComponents::Physical::THREE, IntgPhysType::INTG, IntgPartType::INTG, IntgSpecType::INTG, FieldComponents::Spectral::THREE, FieldType::VECTOR));
      }

      return transform;
   }

   std::vector<ProjectorBranch3D>  backwardScalar(const std::map<FieldComponents::Physical::Id,bool>& req)
   {
      std::vector<ProjectorBranch3D> transform;

      transform.push_back(ProjectorBranch3D(FieldComponents::Spectral::SCALAR, ProjSpecType::PROJ, ProjPartType::PROJ, ProjPhysType::PROJ, FieldComponents::Physical::SCALAR, FieldType::SCALAR));

      return transform;
   }

   std::vector<ProjectorBranch3D>  backwardGradient(const std::map<FieldComponents::Physical::Id,bool>& req)
   {
      std::vector<ProjectorBranch3D> transform;

      if(req.find(FieldComponents::Physical::ONE)->second)
      {
         transform.push_back(ProjectorBranch3D(FieldComponents::Spectral::SCALAR, ProjSpecType::DIFF, ProjPartType::PROJ, ProjPhysType::PROJ, FieldComponents::Physical::ONE, FieldType::GRADIENT));
      }

      if(req.find(FieldComponents::Physical::TWO)->second)
      {
         transform.push_back(ProjectorBranch3D(FieldComponents::Spectral::SCALAR, ProjSpecType::PROJ, ProjPartType::DIFF, ProjPhysType::PROJ, FieldComponents::Physical::TWO, FieldType::GRADIENT));
      }

      if(req.find(FieldComponents::Physical::THREE)->second)
      {
         transform.push_back(ProjectorBranch3D(FieldComponents::Spectral::SCALAR, ProjSpecType::PROJ, ProjPartType::PROJ, ProjPhysType::DIFF, FieldComponents::Physical::THREE, FieldType::GRADIENT));
      }

      return transform;
   }

   std::vector<ProjectorBranch3D>  backwardVector(const std::map<FieldComponents::Physical::Id,bool>& req)
   {
      std::vector<ProjectorBranch3D> transform;

      if(req.find(FieldComponents::Physical::ONE)->second)
      {
         transform.push_back(ProjectorBranch3D(FieldComponents::Spectral::ONE, ProjSpecType::PROJ, ProjPartType::PROJ, ProjPhysType::PROJ, FieldComponents::Physical::ONE, FieldType::VECTOR));
      }

      if(req.find(FieldComponents::Physical::TWO)->second)
      {
         transform.push_back(ProjectorBranch3D(FieldComponents::Spectral::TWO, ProjSpecType::PROJ, ProjPartType::PROJ, ProjPhysType::PROJ, FieldComponents::Physical::TWO, FieldType::VECTOR));
      }

      if(req.find(FieldComponents::Physical::THREE)->second)
      {
         transform.push_back(ProjectorBranch3D(FieldComponents::Spectral::THREE, ProjSpecType::PROJ, ProjPartType::PROJ, ProjPhysType::PROJ, FieldComponents::Physical::THREE, FieldType::VECTOR));
      }

      return transform;
   }

   std::vector<ProjectorBranch3D>  backwardVGradient(FieldComponents::Spectral::Id id, const std::map<FieldComponents::Physical::Id,bool>& req)
   {
      std::vector<ProjectorBranch3D> transform;

      if(req.find(FieldComponents::Physical::ONE)->second)
      {
         transform.push_back(ProjectorBranch3D(id, ProjSpecType::DIFF, ProjPartType::PROJ, ProjPhysType::PROJ, FieldComponents::Physical::ONE, FieldType::GRADIENT));
      }

      if(req.find(FieldComponents::Physical::TWO)->second)
      {
         transform.push_back(ProjectorBranch3D(id, ProjSpecType::PROJ, ProjPartType::DIFF, ProjPhysType::PROJ, FieldComponents::Physical::TWO, FieldType::GRADIENT));
      }

      if(req.find(FieldComponents::Physical::THREE)->second)
      {
         transform.push_back(ProjectorBranch3D(id, ProjSpecType::PROJ, ProjPartType::PROJ, ProjPhysType::DIFF, FieldComponents::Physical::THREE, FieldType::GRADIENT));
      }

      return transform;
   }

   std::vector<ProjectorBranch3D>  backwardCurl(const std::map<FieldComponents::Physical::Id,bool>& req)
   {
      std::vector<ProjectorBranch3D> transform;

      if(req.find(FieldComponents::Physical::ONE)->second)
      {
         transform.push_back(ProjectorBranch3D(FieldComponents::Spectral::TWO, ProjSpecType::PROJ, ProjPartType::PROJ, ProjPhysType::DIFF, FieldComponents::Physical::ONE, FieldType::CURL));

         transform.push_back(ProjectorBranch3D(FieldComponents::Spectral::THREE, ProjSpecType::PROJ, ProjPartType::DIFF, ProjPhysType::PROJ, FieldComponents::Physical::ONE, FieldType::CURL));
      }

      if(req.find(FieldComponents::Physical::TWO)->second)
      {
         transform.push_back(ProjectorBranch3D(FieldComponents::Spectral::ONE, ProjSpecType::PROJ, ProjPartType::PROJ, ProjPhysType::DIFF, FieldComponents::Physical::TWO, FieldType::CURL));

         transform.push_back(ProjectorBranch3D(FieldComponents::Spectral::THREE, ProjSpecType::DIFF, ProjPartType::PROJ, ProjPhysType::PROJ, FieldComponents::Physical::TWO, FieldType::CURL));
      }

      if(req.find(FieldComponents::Physical::THREE)->second)
      {
         transform.push_back(ProjectorBranch3D(FieldComponents::Spectral::ONE, ProjSpecType::PROJ, ProjPartType::DIFF, ProjPhysType::PROJ, FieldComponents::Physical::THREE, FieldType::CURL));

         transform.push_back(ProjectorBranch3D(FieldComponents::Spectral::TWO, ProjSpecType::DIFF, ProjPartType::PROJ, ProjPhysType::PROJ, FieldComponents::Physical::THREE, FieldType::CURL));
      }

      return transform;
   }

   std::vector<ProjectorBranch3D>  backwardDivergence()
   {
      std::vector<ProjectorBranch3D> transform;

      transform.push_back(ProjectorBranch3D(FieldComponents::Spectral::ONE, ProjSpecType::DIFF, ProjPartType::PROJ, ProjPhysType::PROJ, FieldComponents::Physical::SCALAR, FieldType::DIVERGENCE));

      transform.push_back(ProjectorBranch3D(FieldComponents::Spectral::TWO, ProjSpecType::PROJ, ProjPartType::DIFF, ProjPhysType::PROJ, FieldComponents::Physical::SCALAR, FieldType::DIVERGENCE));

      transform.push_back(ProjectorBranch3D(FieldComponents::Spectral::THREE, ProjSpecType::PROJ, ProjPartType::PROJ, ProjPhysType::DIFF, FieldComponents::Physical::SCALAR, FieldType::DIVERGENCE));

      return transform;
   }

}
}
}
